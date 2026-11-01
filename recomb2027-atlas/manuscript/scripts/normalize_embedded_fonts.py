"""Correct CFF programs embedded as FontFile2 by older plotting backends.

Convert only those embedded font programs to TrueType, preserving glyph order,
advance widths, PDF text coordinates and chart geometry. Curves are converted
with at most one font-unit approximation error. Run after regenerating figures.
Requires PyMuPDF and fontTools; no scientific inputs or results are changed.
"""
from io import BytesIO
from fontTools.ttLib import TTFont, newTable
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.pens.cu2quPen import Cu2QuPen
import fitz

def convert(data):
 font=TTFont(BytesIO(data));order=font.getGlyphOrder();gs=font.getGlyphSet();glyphs={}
 for name in order:
  pen=TTGlyphPen(gs);gs[name].draw(Cu2QuPen(pen,1.0,reverse_direction=True));glyphs[name]=pen.glyph()
 del font['CFF ']
 font['glyf']=newTable('glyf');font['glyf'].glyphs=glyphs
 font['loca']=newTable('loca')
 font['maxp']=newTable('maxp');font['maxp'].tableVersion=0x00010000
 for name,value in {'maxZones':1,'maxTwilightPoints':0,'maxStorage':0,'maxFunctionDefs':0,'maxInstructionDefs':0,'maxStackElements':0,'maxSizeOfInstructions':0}.items():setattr(font['maxp'],name,value)
 font['post'].formatType=3.0
 font.sfntVersion='\x00\x01\x00\x00';font.setGlyphOrder(order)
 out=BytesIO();font.save(out);return out.getvalue()

from pathlib import Path

def main():
 root=Path(__file__).resolve().parents[1]
 for path in sorted((root/'figures').glob('*.pdf')):
  d=fitz.open(path);before=[p.get_text() for p in d];changed=0
  for i in range(1,d.xref_length()):
   typ,value=d.xref_get_key(i,'FontFile2')
   if typ=='xref':
    xref=int(value.split()[0]);old=d.xref_stream(xref)
    if old[:4]==b'OTTO':
     data=convert(old);d.update_stream(xref,data)
     d.xref_set_key(xref,'Length1',str(len(data)));changed+=1
  if changed:
   tmp=path.with_suffix('.tmp.pdf');d.save(tmp,garbage=4,deflate=True);d.close()
   check=fitz.open(tmp);assert [p.get_text() for p in check]==before;check.close()
   tmp.replace(path);print(path.name,changed,'font programs corrected')

if __name__=='__main__':main()
