import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
from gudhi import CubicalComplex 
from gudhi import bottleneck_distance

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")

neuron_img = Image.open("mouse_neuron_1_alt.jpg").convert("L")

pixels = np.array(neuron_img)


cc = CubicalComplex(top_dimensional_cells= pixels.max() - pixels)
output = cc.persistence()
h_0 = [[birth, death] for (dim, (birth, death)) in output if dim == 0]

print(bottleneck_distance(h_0, h_0))