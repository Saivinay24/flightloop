"""Index TM-4062 Appendix D (rotated Cp listing) pages by OCR of the header line."""
import subprocess, sys, pathlib
pdf = pathlib.Path(__file__).resolve().parents[1] / "data/tm4062/tm4062.pdf"
work = pdf.parent / "work"; work.mkdir(exist_ok=True)
lo, hi = int(sys.argv[1]), int(sys.argv[2])
for p in range(lo, hi + 1):
    stem = work / f"idx{p}"
    subprocess.run(["pdftoppm", "-f", str(p), "-l", str(p), "-r", "100", "-singlefile", "-png", str(pdf), str(stem)], check=True)
    png = f"{stem}.png"
    subprocess.run(["sips", "-r", "90", png], capture_output=True)
    txt = subprocess.run(["tesseract", png, "-", "--psm", "6"], capture_output=True, text=True).stdout
    hdr = [l.strip() for l in txt.splitlines() if "SURFACE" in l.upper() or "RUN" in l.upper()][:2]
    print(p, " | ".join(hdr), flush=True)
    pathlib.Path(png).unlink()
