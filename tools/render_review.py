"""Render draft PDFs into ignored local contact sheets for visual inspection."""

from pathlib import Path
import subprocess
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
for slug in ("01-projection-sign-rigidity", "02-half-wormholes-wiener-chaos"):
    output = ROOT / "tmp" / slug
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run(["pdftoppm", "-r", "90", "-png", str(ROOT / "papers" / slug / "manuscript.pdf"), str(output / "page")], check=True)
    pages = sorted(output.glob("page-*.png"), key=lambda path: int(path.stem.split("-")[-1]))
    width, height, columns = 320, 475, 3
    sheet = Image.new("RGB", (columns * width, ((len(pages) + columns - 1) // columns) * height), "#dddddd")
    draw = ImageDraw.Draw(sheet)
    for number, path in enumerate(pages):
        with Image.open(path) as page:
            page.thumbnail((width - 12, height - 28))
            x, y = (number % columns) * width, (number // columns) * height
            sheet.paste(page, (x + 6, y + 22))
            draw.text((x + 8, y + 5), f"Page {number + 1}", fill="black")
    target = output / "contact-sheet.jpg"
    sheet.save(target, quality=92)
    print(target)
