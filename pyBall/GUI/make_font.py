# make_font.py
from PIL import Image, ImageDraw, ImageFont
import os

# Configuration
FONT_FILE   = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"  # Ubuntu path
FONT_SIZE   = 32
TILE_W      = 32  # Fixed width per character tile
TILE_H      = 48  # Fixed height per character tile
ASCII_START = 32   # Space ' '
ASCII_END   = 126  # Tilde '~'
NUM_CHARS = ASCII_END - ASCII_START + 1

# 1D horizontal atlas: width = num_chars * tile_width, height = tile_height
IMG_W = NUM_CHARS * TILE_W
IMG_H = TILE_H

def generate_font_atlas():
    """Generate a 1D horizontal font atlas (no JSON needed)."""
    try:
        font = ImageFont.truetype(FONT_FILE, FONT_SIZE)
    except IOError:
        print(f"Font file not found at {FONT_FILE}, using default font.")
        font = ImageFont.load_default( FONT_SIZE )

    img   = Image.new("RGBA", (IMG_W, IMG_H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    for code in range(ASCII_START, ASCII_END + 1):
        char = chr(code)
        x0 = (code - ASCII_START) * TILE_W
        y0 = 0
        draw.text((x0 + 2, 2), char, font=font, fill=(255, 255, 255, 255))

    output_dir = os.path.join(os.path.dirname(__file__), "shaders")
    os.makedirs(output_dir, exist_ok=True)
    img.save(os.path.join(output_dir, "font_atlas.png"))

    # Optional: save a minimal header with only overall dimensions
    with open(os.path.join(output_dir, "font_atlas.json"), "w") as f:
        f.write(f'{{"tile_w":{TILE_W},"tile_h":{TILE_H},"tex_w":{IMG_W},"tex_h":{IMG_H}}}')

if __name__ == "__main__":
    generate_font_atlas()
    
    print("Font atlas 'font_atlas.png' and 'font_atlas.json' created successfully.")

if __name__ == "__main__":
    generate_font_atlas()