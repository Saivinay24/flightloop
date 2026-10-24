"""300 dpi crop of a TM-4062 page region (fractions of width/height), upscaled 2x, for resolving disagreements."""
import sys, subprocess, pathlib
from PIL import Image
ROOT = pathlib.Path(__file__).resolve().parents[1]
work = ROOT / "data/tm4062/work"
page, x0, y0, x1, y1, name = int(sys.argv[1]), *map(float, sys.argv[2:6]), sys.argv[6]
stem = work / f"full{page}"
if not pathlib.Path(f"{stem}.png").exists():
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "300", "-gray", "-singlefile", "-png", str(ROOT / "data/tm4062/tm4062.pdf"), str(stem)], check=True)
im = Image.open(f"{stem}.png"); W, H = im.size
c = im.crop((int(x0 * W), int(y0 * H), int(x1 * W), int(y1 * H)))
c = c.resize((c.width * 2, c.height * 2), Image.LANCZOS)
c.save(work / f"{name}.png"); print(work / f"{name}.png")
