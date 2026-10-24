"""Make 300 dpi crops (x/c column + one value column) of TM-4062 rotated Cp pages for the pass-A manual re-read."""
import sys, subprocess, pathlib
from PIL import Image
ROOT = pathlib.Path(__file__).resolve().parents[1]
work = ROOT / "data/tm4062/work"
page, x0, x1, name = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
stem = work / f"ocr{page}"
if not pathlib.Path(f"{stem}.png").exists():
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "300", "-gray", "-singlefile", "-png", str(ROOT / "data/tm4062/tm4062.pdf"), str(stem)], check=True)
    subprocess.run(["sips", "-r", "90", f"{stem}.png"], capture_output=True)
im = Image.open(f"{stem}.png")
W, H = im.size
xc = im.crop((int(0.10 * W), int(0.25 * H), int(0.17 * W), int(0.86 * H)))
col = im.crop((int(x0 * W), int(0.25 * H), int(x1 * W), int(0.86 * H)))
out = Image.new("L", (xc.width + col.width, xc.height), 255)
out.paste(xc, (0, 0)); out.paste(col, (xc.width, 0))
out.save(work / f"{name}.png")
print(W, H, work / f"{name}.png")
