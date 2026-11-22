# Upstream Direct-Prior MIL Reference Code (`references/code/`)

To guarantee self-contained in-repository reproducibility without relying on out-of-tree `~/.cache/crlm-m4-priors/` clones, the exact upstream source trees for the two direct-prior multiple-instance learning (MIL) methods referenced in `docs/CRLM_Research_Plan.md` (`[M01]`, `[M02]`) and `src/crlm/prior_mil.py` are archived here at fixed commits (with `.git/` stripped to prevent nested Git repository conflicts).

## 1. Archived Upstream Repositories

| Directory | Reference ID | Upstream URL | Fixed Commit | Commit Date | License | Files | Size |
|---|---|---|---|---|---|---:|---:|
| `references/code/hier-mil/` | `[M02]` Do & Lähdesmäki (2025) | `https://github.com/minhchaudo/hier-mil` | `2fb155840394961595755bede506c4b54c0942ae` | `2025-04-10T10:33:37+03:00` | No explicit LICENSE file in upstream repository (public academic reference code archived read-only for audit/reproduction) | 18 | 29,685,158 B |
| `references/code/mixmil/` | `[M01]` Engelmann, Palma & Casale (2024/2025) | `https://github.com/AIH-SGML/mixmil` | `16d92349ccf19b73c1dccd167934fcc3202b6be4` | `2025-11-10T18:13:22+01:00` | `Apache-2.0 (see references/code/mixmil/LICENSE)` | 17 | 76,172 B |

## 2. Key File SHA256 Checksums

### `hier-mil` (`2fb155840394`)
- `references/code/hier-mil/README.md`: `a5abb30db804a0fbb19310c4c9a969089c65d8eec5841a0cf7cf88306bd26de4`
- `references/code/hier-mil/model.py`: `fefdb233563b56de3030e65d0b8ef203bb2dca19c342d9bc87defb6b2f460e2c`
- `references/code/hier-mil/requirements.txt`: `8555b082e32a7c55b0946bfd5f003e810295d18a827138d62da948910d06eed9`
- `references/code/hier-mil/run.py`: `f09f0e979c49bc8732082c12ca8461caa4d40641fde696007fb8ff04dd7362a0`
- `references/code/hier-mil/train.py`: `ac2c72ad7207253b394cafb23f122e621a506cab9d67ea2d0663798fa9fe1237`
- `references/code/hier-mil/utils.py`: `d1323b1927a79788430af448e904d99054033b7a33d7840a46f393c213d1696e`

### `mixmil` (`16d92349ccf1`)
- `references/code/mixmil/LICENSE`: `c71d239df91726fc519c6eb72d318ec65820627232b2f796219e87dcf35d0ab4`
- `references/code/mixmil/README.md`: `24317cf23f582ef52b093fda148a2aa68a3687e7ec93e80b66ea2c92377da5b5`
- `references/code/mixmil/mixmil/__init__.py`: `7d77f273e2176317d43c842d1140bc1d2b76d03c6c6b6b533e2f604cefb81c2e`
- `references/code/mixmil/mixmil/data.py`: `2d8ee2ad3a08c13053f14cee311fcbaf687dd2eca641b7ad15b39fca73f5bfb8`
- `references/code/mixmil/mixmil/model.py`: `be10e3a865c21c420d9d8a97b4dd9e1c175e9de789f5556d227daa8152062310`
- `references/code/mixmil/mixmil/posterior.py`: `64c352ae4127526f1ff2952305484f3b3cabce05ade29e7af61e146775f9ce89`
- `references/code/mixmil/mixmil/utils.py`: `a96ba88d4f54040c25334ca3d7b4c6889e46c1d51073b9561bb632e1d56c6203`
- `references/code/mixmil/pyproject.toml`: `99ff918d77654b935c6adde6f95a12d81d0562ca3e0416eec69922be6288f446`

## 3. Usage & Audit Notes

- Full per-file SHA256 inventory (all 35 archived files across both trees) is recorded in `references/code/PROVENANCE.json`.
- Project benchmark scripts (`src/crlm/prior_mil.py`, `scripts/benchmark_prior_mil.py`) align against these exact commits on CPU/XPU without requiring `torch_scatter` CUDA extensions (using pure PyTorch index operations equivalent to `mixmil/model.py` and `hier-mil/model.py`).
