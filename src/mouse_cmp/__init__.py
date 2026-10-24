"""mouse_cmp — mouse DMG model vs human DMG atlas comparison, and treatment (FUS / RT / FUS+RT / ICI) effects.

Shares thread caps and StepLogger with the main pipeline: scripts put both `mouse_comparison/src` and the repo `src/`
on sys.path and `import dmg_atlas.utils` before numpy/scanpy.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # this repo (dmg_mh_states)
# The human atlas was a parent directory when this sub-project lived inside it; it is now a sibling repository.
# Resolve it explicitly: $DMG_ATLAS_ROOT wins, else a sibling 'dmg_atlas', else the old parent layout.
import os as _os
_env = _os.environ.get("DMG_ATLAS_ROOT")
if _env:
    ATLAS_ROOT = Path(_env)
elif (ROOT.parent / "dmg_atlas" / "data" / "processed").is_dir():
    ATLAS_ROOT = ROOT.parent / "dmg_atlas"
else:
    ATLAS_ROOT = ROOT.parent
RAW = ROOT / "data" / "raw"                         # mouse deposits (Seurat RDS as received)
INTERIM = ROOT / "data" / "interim"                 # converted / QC'd mouse h5ad
PROCESSED = ROOT / "data" / "processed"             # ortholog-mapped, annotated mouse objects
HUMAN = ROOT / "data" / "human"                     # read-only symlinks into the human atlas outputs
EXTERNAL = ROOT / "data" / "external"               # ortholog tables, mouse reference gene sets
RESULTS = ROOT / "results"
CONFIGS = ROOT / "configs"
LOGS = ROOT / "logs"
