"""Static HTML export of the dashboard: the self-contained site/ (index.html + assets + data/bundle.js).
Works from file:// and under any subpath (relative paths only)."""
import runpy, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
runpy.run_path(str(ROOT / "scripts/build_site.py"), run_name="__main__")
print(f"static export: {ROOT / 'site/index.html'}")
