import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
import gudhi as ghd 
import MP_image_processing as mp

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")



pixels = np.zeros((5,5))
ring = mp.nth_ring(pixels, 2, 2, 1)
for k,l in ring:
    pixels[k,l] = 1

st = ghd.SimplexTree()
mp.insert_simplices(st, pixels, mp.top)
print(st.num_simplices())

for l,f in st.get_simplices(): 
    print("Simplex found:", l, "with filtration value", f)


st.compute_persistence()
pers_list = st.persistence() 

ghd.plot_persistence_diagram(st.persistence())
plt.show()