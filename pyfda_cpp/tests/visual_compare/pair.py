import sys, os
from PIL import Image, ImageDraw
sys.path.insert(0, os.path.dirname(__file__))
from designs import D
root = sys.argv[1]; out = os.path.join(root, 'pairs'); os.makedirs(out, exist_ok=True)
PAIRS = [('H', 'plot0', 'tab0'), ('phi', 'plot1', 'tab1'), ('tau', 'plot2', 'tab2'), ('pz', 'plot3', 'tab3'),
         ('3d', 'plot5', 'tab4'), ('yn', 'plot4', 'tab6'), ('coeffs', 'input1', 'tab7'), ('fix', 'input4', 'tab8'),
         ('info', 'input3', 'tab9'), ('datafilt', 'plot6', 'tab10')]
for name, _ in D:
    for key, p, c in PAIRS:
        if key == 'datafilt' and name != 'd1_ellip_lp': continue
        a = Image.open(f'{root}/py/{name}/{p}.png').convert('RGB'); b = Image.open(f'{root}/cpp/{name}/{c}.png').convert('RGB')
        h = 820
        a = a.resize((int(a.width * h / a.height), h)); b = b.resize((int(b.width * h / b.height), h))
        im = Image.new('RGB', (a.width + b.width + 10, h + 28), 'white')
        im.paste(a, (0, 28)); im.paste(b, (a.width + 10, 28))
        d = ImageDraw.Draw(im); d.text((6, 8), 'pyfda (Python)', fill='black'); d.text((a.width + 16, 8), 'pyfda C++', fill='black')
        im.save(f'{out}/{name}__{key}.png')
