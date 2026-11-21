"""Count donor/organ support in the original healthy-organ atlas annotation.

Uses only obs obtained by probe_normal_atlas_metadata.py. Does not reclassify
cells, read expression, or estimate biological effects. Caecum is kept separate
from transverse/sigmoid colon and small intestine.
"""
from collections import Counter
from pathlib import Path
import csv
import json
import gzip

ROOT = Path(__file__).resolve().parent
obs_path = ROOT/"geo/E-MTAB-11536_obs_selected.json"
obs = json.loads(obs_path.read_text() if obs_path.exists() else
                 gzip.decompress(obs_path.with_suffix(".json.gz").read_bytes()).decode())
counts = Counter(zip(obs["Donor"], obs["Organ"]))
rows = [dict(donor=d, organ=o, n_immune_cells=n) for (d,o),n in sorted(counts.items())]
with (ROOT/"E-MTAB-11536_donor_organ_cell_counts.csv").open("w", newline="") as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
lineage_counts = Counter()
for donor, organ, celltype in zip(obs["Donor"],obs["Organ"],obs["Majority_voting_CellTypist_high"]):
    if organ not in {"LIV", "CAE", "TCL", "SCL"}: continue
    lineage = ("Myeloid" if celltype in {"Monocytes","Macrophages","DC","pDC"} else
               "T_NK" if celltype in {"T cells","ILC"} else
               "B_plasma" if celltype in {"B cells","Plasma cells","B-cell lineage"} else "Other")
    lineage_counts[(donor,organ,lineage)] += 1
rows = [dict(donor=d,organ=o,lineage=l,n_cells=n) for (d,o,l),n in sorted(lineage_counts.items())]
with (ROOT/"E-MTAB-11536_liver_colon_lineage_support.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
summary={
    "n_immune_cells":len(obs["Donor"]), "n_donors":len(set(obs["Donor"])),
    "liver_donors":sorted({d for d,o in counts if o=="LIV"}),
    "transverse_sigmoid_donors":sorted({d for d,o in counts if o in {"TCL","SCL"}}),
    "caecum_donors":sorted({d for d,o in counts if o=="CAE"}),
    "annotation_source":"original global.h5ad Majority_voting_CellTypist_high; not reannotation",
    "limitation":"Not all12donors havecolon/liver. Caecum, transversecolon, sigmoidcolon and smallintestine are distinct sites; sparsemyeloidcolon prevents pairedmyeloidnormal inference at50cells/site.",
}
(ROOT/"E-MTAB-11536_support_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary))
