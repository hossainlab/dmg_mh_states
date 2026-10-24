"""mouse_cmp — mouse DMG model vs human DMG atlas comparison, and treatment (FUS / RT / FUS+RT / ICI) effects.

Shares thread caps and StepLogger with the main pipeline: scripts put both `mouse_comparison/src` and the repo `src/`
on sys.path and `import dmg_atlas.utils` before numpy/scanpy.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # mouse_comparison/
ATLAS_ROOT = ROOT.parent                            # dmg_atlas repo root
RAW = ROOT / "data" / "raw"                         # mouse deposits (Seurat RDS as received)
INTERIM = ROOT / "data" / "interim"                 # converted / QC'd mouse h5ad
PROCESSED = ROOT / "data" / "processed"             # ortholog-mapped, annotated mouse objects
HUMAN = ROOT / "data" / "human"                     # read-only symlinks into the human atlas outputs
EXTERNAL = ROOT / "data" / "external"               # ortholog tables, mouse reference gene sets
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs"
LOGS = ROOT / "logs"
