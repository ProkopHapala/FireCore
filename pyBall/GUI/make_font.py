# make_font.py
from PIL import Image, ImageDraw, ImageFont
import json
import os

# Configuration
FONT_FILE = "C:/Windows/Fonts/consola.ttf" # Or a path to any other monospace font
FONT_SIZE = 32
CHARS = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
IMG_WIDTH = 512
PADDING = 2

def generate_font_atlas():
    """Generates a font atlas image and a JSON file with character metrics."""
    try:
        font = ImageFont.truetype(FONT_FILE, FONT_SIZE)
    except IOError:
        print(f"Font file not found at {FONT_FILE}, using default font.")
        font = ImageFont.load_default()

    img = Image.new("RGBA", (IMG_WIDTH, 1024), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    char_data = {}
    x, y, max_h = PADDING, PADDING, 0

    for char in CHARS:
        # Use getbbox for modern Pillow versions
        bbox = draw.textbbox((0,0), char, font=font)
        char_width = bbox[2] - bbox[0]
        char_height = bbox[3] - bbox[1]

        if x + char_width + PADDING > IMG_WIDTH:
            x = PADDING
            y += max_h + PADDING
            max_h = 0
        
        if char_height > max_h:
            max_h = char_height

        draw.text((x, y), char, font=font, fill=(255, 255, 255, 255))
        
        char_data[char] = {
            'x': x, 'y': y,
            'w': char_width, 'h': char_height,
            'norm_x': x / IMG_WIDTH, 'norm_y': y / 1024,
            'norm_w': char_width / IMG_WIDTH, 'norm_h': char_height / 1024,
        }
        x += char_width + PADDING

    # Crop image to actual height
    final_img = img.crop((0, 0, IMG_WIDTH, y + max_h + PADDING))
    output_dir = os.path.join(os.path.dirname(__file__), "shaders")
    os.makedirs(output_dir, exist_ok=True)
    final_img.save(os.path.join(output_dir, "font_atlas.png"))

    # Save metrics
    with open(os.path.join(output_dir, "font_atlas.json"), "w") as f:
        json.dump({'chars': char_data, 'tex_w': final_img.width, 'tex_h': final_img.height}, f)
    
    print("Font atlas 'font_atlas.png' and 'font_atlas.json' created successfully.")

if __name__ == "__main__":
    generate_font_atlas()