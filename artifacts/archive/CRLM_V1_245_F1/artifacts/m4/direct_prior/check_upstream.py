"""Numerically compare adapted pooling to the pinned author implementations."""
import ast
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"src"))
from crlm.prior_mil import BagBatch, HierarchicalMIL, MixedMIL, official_posterior


def verify():
    cache=Path.home()/".cache/crlm-m4-priors"
    revisions={k:subprocess.check_output(["git","-C",str(cache/k),"rev-parse","HEAD"],text=True).strip()
               for k in ("mixmil","hier-mil")}
    spec=importlib.util.spec_from_file_location("author_hier",cache/"hier-mil/model.py")
    source=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    torch.manual_seed(123)
    raw=np.random.default_rng(23).normal(size=(8,8,4)).astype("f")
    labels=[np.array(["a"]*4+["b"]*4) for _ in raw]
    batch=BagBatch(raw,labels)
    differences={}
    for attention in (False,True):
        ref=source.Model(4,attn1=attention,attn2=True,dropout=0,n_layers_lin=1,n_hid=32)
        model=HierarchicalMIL(4,attention)
        model.transform.load_state_dict(ref.lin[0].state_dict())
        model.type_attention.load_state_dict(ref.w_ct[0].state_dict())
        model.head.load_state_dict(ref.lin_out.state_dict())
        if attention:
            model.cell_attention.load_state_dict(ref.w_c[0].state_dict())
        expected=ref(batch.x,batch.group,16,2).flatten()
        differences["HA" if attention else "CTA"]=(expected-model(batch)[0]).abs().max().item()

    # The published tensor kernel only depends on torch; use that exact AST
    # without importing the unrelated torch_scatter path for ragged batches.
    mix_path=cache/"mixmil/mixmil/model.py"
    module=ast.parse(mix_path.read_text())
    cls=next(c for c in module.body if isinstance(c,ast.ClassDef) and c.name=="MixMIL")
    method=next(f for f in cls.body if isinstance(f,ast.FunctionDef) and f.name=="_calc_bag_emb_effect_tensor")
    symbols={"torch":torch}
    exec(compile(ast.Module(body=[method],type_ignores=[]),str(mix_path),"exec"),symbols)
    model=MixedMIL(raw,np.array([0,1]*4),official_posterior(cache/"mixmil"))
    beta=torch.randn(8,8)
    pooled,_,_,_=model.raw(batch,beta)
    u,z=beta.split(4)
    eta=z/z.square().mean(0).sqrt()
    expected=symbols[method.name](None,u[:,None,:],eta[:,None,:],torch.tensor(raw)).squeeze(1)
    differences["MixMIL"]=(expected-pooled).abs().max().item()
    if max(differences.values())>2e-6:
        raise ValueError(differences)
    dest=Path(__file__).parent/"upstream_equivalence.json"
    dest.write_text(json.dumps({"revisions":revisions,"max_absolute_difference":differences,
        "scope":"all subtypes observed for HA/CTA; same beta and raw pre-scaling bag effects for MixMIL"},indent=2)+"\n")


if __name__=="__main__":
    verify()
