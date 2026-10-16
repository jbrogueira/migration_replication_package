# tests/test_paper_exhibits.py
# The paper's table and figure files are built from existing outputs (run_paper.py);
# skipped when those outputs are absent (fresh checkout before the full run).
import filecmp
import os

import pytest

from analysis import build_paper_figures, build_paper_tables

PKG = build_paper_tables.PKG
OUT = os.path.join(PKG, "output")
NEEDED = [os.path.join(OUT, "tables", "tex"), os.path.join(OUT, "figures"),
          *(os.path.join(OUT, "variants", v) for v in ("eps20", "gamma005", "eu_born", "sigma_theta_grid"))]
pytestmark = pytest.mark.skipif(not all(os.path.isdir(p) for p in NEEDED),
                                reason="outputs of run_paper.py not found")


def test_build_paper_tables(tmp_path):
    written = build_paper_tables.build(str(tmp_path))
    assert len(written) >= len(build_paper_tables.tables())
    for name, text in written.items():
        assert text.lstrip().startswith(r"\begin{tabular}"), name
        assert text.endswith("\\end{tabular}%\n"), name
        assert (tmp_path / f"{name}.tex").read_text(encoding="utf-8") == text
    # inside the paper repository, the draft's table files are the generated ones
    paper = os.path.normpath(os.path.join(PKG, "..", "..", "text", "tables"))
    if os.path.isdir(paper):
        stale = [n for n in written
                 if not filecmp.cmp(tmp_path / f"{n}.tex", os.path.join(paper, f"{n}.tex"), shallow=False)]
        assert not stale, stale


def test_build_paper_figures(tmp_path):
    src = build_paper_figures.build(str(tmp_path))
    assert set(build_paper_figures.STATIC) <= set(src)
    for name, path in src.items():
        assert filecmp.cmp(path, tmp_path / name, shallow=False)
