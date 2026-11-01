"""Reletter the recorded donor map at its final 6.5-inch publication width.

This preserves the original points, region boxes, axes and coordinates. Only
text is replaced. No tissue data are reconstructed from the rendered image.
Requires PyMuPDF and Nimbus Sans (URW Base35).
"""
from pathlib import Path
import fitz

ROOT = Path(__file__).resolve().parents[1]
FONT_DIR = Path('/usr/share/fonts/opentype/urw-base35')

def main():
    original = fitz.open(ROOT / 'figure_data/biological_regions_original.pdf')
    source = original[0]
    spans = [(line['dir'], span) for block in source.get_text('dict')['blocks']
             for line in block.get('lines', []) for span in line['spans']]
    for _, span in spans:
        source.add_redact_annot(fitz.Rect(span['bbox']), fill=None, cross_out=False)
    source.apply_redactions(images=0, graphics=0, text=0)
    factor = 468 / source.rect.width
    offset = 10
    out = fitz.open()
    page = out.new_page(width=468, height=source.rect.height * factor + 20)
    page.show_pdf_page(fitz.Rect(0, offset, 468, source.rect.height * factor + offset), original, 0)
    fonts = {}
    for key, name in [('regular', 'NimbusSans-Regular'), ('bold', 'NimbusSans-Bold')]:
        fontpath = FONT_DIR / (name + '.otf')
        page.insert_font(fontname=key, fontfile=str(fontpath))
        fonts[key] = fitz.Font(fontfile=str(fontpath))
    for direction, span in spans:
        text = span['text']
        if text == '1e3':
            continue
        if text in ('x (µm)', 'y (µm)'):
            text = text.replace('µm', 'mm')
        rect = fitz.Rect(span['bbox'])
        cx = (rect.x0 + rect.x1) * factor / 2
        cy = (rect.y0 + rect.y1) * factor / 2 + offset
        is_title = 'carcinoma donor' in text
        is_region = 'Bold' in span['font'] and not is_title
        key = 'bold' if is_title or is_region else 'regular'
        size = 10.5 if is_title else 10.0
        if is_title:
            text = text.replace(' carcinoma', '')
        elif text == '800 µm context':
            text = 'Context'
        elif text == '200 µm core':
            text = 'Core'
        font = fonts[key]
        width = font.text_length(text, fontsize=size)
        rgb = tuple(((span['color'] >> shift) & 255) / 255 for shift in (16, 8, 0))
        if direction == (0.0, -1.0):
            # Keep a small clear left margin for the rotated y-axis title.
            x = max(8.4, cx) + size * (font.ascender + font.descender) / 2
            point = (x, cy + width / 2)
            page.insert_text(point, text, fontsize=size, fontname=key, color=rgb, rotate=90)
        else:
            x = cx - width / 2
            if span['text'] in ('800 µm context', '200 µm core'):
                x = rect.x0 * factor
            # Left/right y ticks align against the same axes as the original.
            elif text.isdigit() and not is_region and rect.x1 < 35:
                x = rect.x1 * factor - width
            elif text.isdigit() and not is_region and 310 < rect.x0 < 323:
                x = rect.x1 * factor - width
            x = min(max(x, 1), 467 - width)
            y = cy + size * (font.ascender + font.descender) / 2
            page.insert_text((x, y), text, fontsize=size, fontname=key, color=rgb)
    target = ROOT / 'figures/biological_regions.pdf'
    out.set_metadata({'title':'Predefined breast cancer donor regions', 'author':'Arnav Peram'})
    out.save(target, garbage=4, deflate=True)
    page.get_pixmap(matrix=fitz.Matrix(2,2)).save(ROOT / 'figures/biological_regions.png')
    (ROOT / 'figures/biological_regions.svg').write_text(page.get_svg_image(text_as_path=False))
    sizes = [s['size'] for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    assert min(sizes) >= 10, min(sizes)
    print(f'{target}: {len(sizes)} labels, minimum {min(sizes):g} PDF points')

if __name__ == '__main__':
    main()
