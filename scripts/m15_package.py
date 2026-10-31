"""Build the collaborator package: figures, report and index in one self-contained folder.

Why this exists: REPORT.md links its figures as `../results/figures/<name>.png`, which is correct inside the repo
and broken everywhere else. Zipping the report and the figures into a flat package left all thirteen image links
dangling, so a recipient opened the report to thirteen broken images. The packaged copy of the report therefore
has its links rewritten to the package's own layout, and the PNGs travel with it.

Outputs: dist/dmg_mouse_human_report.zip
Run: python scripts/m15_package.py
"""
from __future__ import annotations

import re
import shutil
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS, DOCS, DIST = ROOT / "results", ROOT / "docs", ROOT / "dist"
PKG = DIST / "pkg"
ZIP = DIST / "dmg_mouse_human_report.zip"


def main() -> int:
    if PKG.exists():
        shutil.rmtree(PKG)
    (PKG / "figures").mkdir(parents=True)
    (PKG / "tables").mkdir(parents=True)

    for src in sorted(RESULTS.glob("figures/*.pdf")) + sorted(RESULTS.glob("figures_denovo/*.pdf")) + \
               sorted(RESULTS.glob("figures/*.png")) + sorted(RESULTS.glob("figures_denovo/*.png")):
        shutil.copy2(src, PKG / "figures" / src.name)
    shutil.copy2(RESULTS / "tables/m13_figure_index.csv", PKG / "tables" / "m13_figure_index.csv")

    # Rewrite the report's repo-relative figure links to the package layout.
    text = (DOCS / "REPORT.md").read_text()
    text, n = re.subn(r"\.\./results/figures(?:_denovo)?/", "figures/", text)
    (PKG / "REPORT.md").write_text(text)

    missing = [m for m in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text) if not (PKG / m).exists()]
    if missing:
        print("BROKEN LINKS IN PACKAGE:", missing)
        return 1

    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(PKG.rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(DIST))
    print(f"{ZIP.name}: {len(list(PKG.glob('figures/*.pdf')))} PDFs, "
          f"{len(list(PKG.glob('figures/*.png')))} PNGs, {n} figure links rewritten, 0 broken; "
          f"{ZIP.stat().st_size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
