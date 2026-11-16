import numpy as np, sys
from PIL import Image
im=np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
dark=(im.sum(2)<200)
H,W=dark.shape
cols=dark.sum(0); rows=dark.sum(1)
print("long vertical dark lines (x, count):",[(x,c) for x,c in enumerate(cols) if c>0.5*H])
print("long horizontal dark lines (y, count):",[(y,c) for y,c in enumerate(rows) if c>0.3*W])
