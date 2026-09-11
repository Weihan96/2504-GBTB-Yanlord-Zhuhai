"""Render existing native SVG evidence while independent schema checking runs."""
from pathlib import Path
import subprocess,sys
from PIL import Image,ImageOps,ImageDraw
OUT=Path(__file__).resolve().parent
with_single='--with-single' in sys.argv
board=Image.new('RGB',(1800,1500 if with_single else 1000),'white');draw=ImageDraw.Draw(board)
for i,view in enumerate(['PLAN','FRONT','SIDE']):
 svg=OUT/f'STREET-H-SUPPORT-SCENE-{view}.svg';png=svg.with_suffix('.png')
 if not with_single or not png.exists():subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1','--export-width=1300',f'--export-filename={png}'],check=True,capture_output=True)
 im=ImageOps.contain(Image.open(png).convert('RGB'),(570,900))
 if with_single:im=ImageOps.contain(im,(570,800))
 board.paste(im,(i*600+(600-im.width)//2,650 if with_single else 60));draw.text((i*600+25,615 if with_single else 25),f'STREET-H SUPPORT | {view} | SCENE PENDING',fill='black')
 if with_single:
  single=ImageOps.contain(Image.open(OUT/f'STREET-H-SUPPORT-SINGLE-{view}.png').convert('RGB'),(560,500))
  board.paste(single,(i*600+(600-single.width)//2,60));draw.text((i*600+25,25),f'STREET-H SUPPORT | {view} | APPROVED SINGLE',fill='black')
board.save(OUT/('STREET-H-SUPPORT-contact-sheet.png' if with_single else 'STREET-H-SUPPORT-scene-contact-sheet.png'))
