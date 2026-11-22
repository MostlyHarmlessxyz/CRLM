"""Patient-weighted adapters for hierarchical MIL and full-posterior MixMIL."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.distributions import LowRankMultivariateNormal, kl_divergence
from torch.nn import functional as F

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MIXMIL_DIR = ROOT / "references/code/mixmil"
DEFAULT_HIERMIL_DIR = ROOT / "references/code/hier-mil"


def resolve_prior_repo(name: str, source: Path | str | None = None) -> Path:
    """Resolve an archived upstream prior repository inside references/code/."""
    defaults = {
        "mixmil": DEFAULT_MIXMIL_DIR,
        "hier-mil": DEFAULT_HIERMIL_DIR,
    }
    if name not in defaults:
        raise ValueError(f"Unknown prior repository {name!r}; expected one of {sorted(defaults)}")
    repo = defaults[name] if source is None else Path(source)
    if not repo.is_dir():
        raise FileNotFoundError(f"Missing {name} source directory: {repo}")
    return repo


class BagBatch:
    """Flatten cells into observed subtypes within individual specimens."""

    def __init__(self, matrices, types):
        self.x = torch.as_tensor(np.concatenate(matrices), dtype=torch.float32)
        self.n = len(matrices)
        self.owner = torch.as_tensor(np.repeat(np.arange(self.n), [len(x) for x in matrices]))
        group_index, group_owner = [], []
        for i, labels in enumerate(types):
            names, inv = np.unique(labels, return_inverse=True)
            group_index.extend(inv + len(group_owner))
            group_owner.extend([i] * len(names))
        self.group = torch.as_tensor(group_index)
        self.group_owner = torch.as_tensor(group_owner)


def group_sum(x, index, n):
    return x.new_zeros((n, x.shape[1])).index_add(0, index, x)


def group_softmax(x, index, n):
    indices = index[:, None].expand_as(x)
    maximum = x.new_full((n, x.shape[1]), -torch.inf)
    maximum.scatter_reduce_(0, indices, x.detach(), reduce="amax")
    val = (x - maximum[index]).exp()
    return val / group_sum(val, index, n)[index]


class HierarchicalMIL(nn.Module):
    """Do–Lähdesmäki one-layer HA/CTA with explicit observed-type masks and optional dropout."""

    def __init__(
        self,
        width: int,
        cell_attention: bool,
        *,
        hidden_dim: int = 32,
        dropout: float = 0.0,
    ):
        super().__init__()
        if hidden_dim < 1 or not (0.0 <= float(dropout) < 1.0):
            raise ValueError("Require hidden_dim >= 1 and 0.0 <= dropout < 1.0")
        self.hidden_dim = int(hidden_dim)
        self.dropout_p = float(dropout)
        self.transform = nn.Linear(width, self.hidden_dim)
        self.drop = nn.Dropout(self.dropout_p)
        self.cell_attention = nn.Linear(self.hidden_dim, 1) if cell_attention else None
        self.type_attention = nn.Linear(self.hidden_dim, 1)
        self.head = nn.Linear(self.hidden_dim, 1)

    def forward(self, bags):
        h = self.drop(self.transform(bags.x).relu())
        n_groups = len(bags.group_owner)
        if self.cell_attention is not None:
            a = group_softmax(self.drop(self.cell_attention(h)), bags.group, n_groups)
        else:
            a = 1 / torch.bincount(bags.group)[bags.group, None].to(h.dtype)
        ht = group_sum(h * a, bags.group, n_groups)
        b = group_softmax(self.drop(self.type_attention(ht)), bags.group_owner, bags.n)
        logits = self.head(group_sum(ht * b, bags.group_owner, bags.n)).flatten()
        weights = (a * b[bags.group]).flatten()
        signed = weights * (h @ self.head.weight.T).flatten()
        return logits, weights, signed, self.head.bias[0]


def official_posterior(source: Path | str | None = None):
    repo = resolve_prior_repo("mixmil", source)
    path = repo / "mixmil/posterior.py"
    if not path.is_file():
        raise FileNotFoundError(f"Missing MixMIL posterior module: {path}")
    spec = importlib.util.spec_from_file_location("mixmil_upstream_posterior", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.GaussianVariationalPosterior


def official_hiermil_model(source: Path | str | None = None):
    repo = resolve_prior_repo("hier-mil", source)
    path = repo / "model.py"
    if not path.is_file():
        raise FileNotFoundError(f"Missing hier-mil model module: {path}")
    spec = importlib.util.spec_from_file_location("hiermil_upstream_model", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Model


def mean_glmm_init(matrices, labels, *, legacy_even_slice: bool = False):
    from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

    x = np.stack([a.mean(axis=0) for a in matrices]).astype(float)
    x -= x.mean(axis=0)
    x /= np.sqrt((x * x).sum() / (len(x) - 1))
    model = BinomialBayesMixedGLM(labels, np.ones((len(x), 1)), x, np.zeros(x.shape[1], dtype=int))
    fit = model.fit_vb()
    # Upstream mixmil/utils.py:_get_single_binomial_init_params sliced [::2] to undo
    # Xm.repeat(2, axis=0) under n_trials=2, leaving a latent half-sample bug when
    # n_trials=1 (no row duplication). Use all training bags unless reproducing the
    # archived pre-remediation baseline via legacy_even_slice=True.
    proj = x @ fit.vc_mean
    if legacy_even_slice:
        proj = proj[::2]
    denom = max(float(np.sqrt(np.mean(fit.vc_mean**2))), 1e-12)
    ratio = float(np.std(proj)) / denom
    return tuple(
        torch.as_tensor(v, dtype=torch.float32)
        for v in (ratio * fit.vc_mean[:, None], ratio * fit.vc_sd[:, None], fit.fe_mean[0])
    )


class MixedMIL(nn.Module):
    """Full-covariance Gaussian variational MixMIL adapter with optional 1/sqrt(Q) scaling."""

    def __init__(
        self,
        matrices,
        labels,
        posterior_class=None,
        *,
        scale_by_sqrt_q: bool = True,
        legacy_even_slice_init: bool = False,
    ):
        super().__init__()
        if posterior_class is None:
            posterior_class = official_posterior()
        mean, sd, intercept = mean_glmm_init(
            matrices, labels, legacy_even_slice=legacy_even_slice_init
        )
        self.q = len(mean)
        self.scale_by_sqrt_q = bool(scale_by_sqrt_q)
        feature_scale = 1.0 / math.sqrt(self.q) if self.scale_by_sqrt_q else 1.0
        self.posterior = posterior_class(2 * self.q, 1, False, init_params=(mean, sd))
        self.log_sigma_u = nn.Parameter(torch.tensor(0.5 * np.log(0.5), dtype=torch.float32))
        self.log_sigma_z = nn.Parameter(
            0.5 * torch.log(torch.mean(mean**2 + sd**2).clamp_min(1e-12))
        )
        self.intercept = nn.Parameter(intercept)
        self.register_buffer("feature_scale", torch.tensor(feature_scale, dtype=torch.float32))
        self.register_buffer("training_mean", torch.zeros(1))
        self.register_buffer("training_sd", torch.ones(1))

    def raw(self, bags, beta):
        u, z = beta.split(self.q, dim=0)
        b = torch.sqrt(z.square().mean(dim=0).clamp_min(1e-12))
        eta = z / b
        x_scaled = bags.x * self.feature_scale
        a = group_softmax(x_scaled @ u, bags.owner, bags.n)
        contribution = a * (x_scaled @ eta)
        return group_sum(contribution, bags.owner, bags.n), a, contribution, b

    def variational(self, bags):
        posterior = self.posterior.distribution
        beta = posterior.rsample((8,)).squeeze(1).T
        pooled, _, _, scale = self.raw(bags, beta)
        logits = self.intercept + scale * (pooled - pooled.mean(0)) / pooled.std(0).clamp_min(1e-8)
        diagonal = (
            torch.cat((self.log_sigma_u.expand(self.q), self.log_sigma_z.expand(self.q)))
            .exp()
            .square()[None, :]
        )
        prior = LowRankMultivariateNormal(
            torch.zeros_like(diagonal), torch.zeros_like(diagonal[..., None]), diagonal
        )
        return logits, kl_divergence(posterior, prior).sum()

    @torch.no_grad()
    def bind_training_scale(self, bags):
        pooled, _, _, _ = self.raw(bags, self.posterior.q_mu)
        self.training_mean.copy_(pooled.mean(0))
        self.training_sd.copy_(pooled.std(0).clamp_min(1e-8))

    def forward(self, bags):
        pooled, attention, cell, b = self.raw(bags, self.posterior.q_mu)
        signed = b * cell / self.training_sd
        bias = self.intercept - b * self.training_mean / self.training_sd
        logits = bias + b * pooled / self.training_sd
        return logits.flatten(), attention.flatten(), signed.flatten(), bias[0]


def train_fixed(
    model,
    bags,
    design,
    n_donors,
    steps,
    *,
    lr: float = 0.001,
    weight_decay: float | None = None,
    lambda_pair: float = 0.1,
):
    weights = torch.tensor(design.bce_weights, dtype=torch.float32)
    labels = torch.tensor(design.labels, dtype=torch.float32)
    contrasts = torch.tensor(design.pair_matrix, dtype=torch.float32)
    mixed = isinstance(model, MixedMIL)
    wd = (0.0 if mixed else 0.001) if weight_decay is None else float(weight_decay)
    optim = torch.optim.Adam(model.parameters(), lr=float(lr), weight_decay=wd)
    model.train()
    loss = torch.tensor(0.0)
    for _ in range(steps):
        if mixed:
            logits, kl = model.variational(bags)
        else:
            logits = model(bags)[0][:, None]
            kl = 0.0
        bce = F.binary_cross_entropy_with_logits(
            logits, labels[:, None].expand_as(logits), reduction="none"
        )
        pair_term = (
            float(lambda_pair) * F.softplus(-(contrasts @ logits)).mean()
            if len(contrasts) and lambda_pair > 0
            else 0.0
        )
        loss = (weights @ bce).mean() + pair_term + kl / n_donors
        if not torch.isfinite(loss):
            raise FloatingPointError("Non-finite training loss")
        optim.zero_grad()
        loss.backward()
        optim.step()
    model.eval()
    if mixed:
        model.bind_training_scale(bags)
    return loss.item()
