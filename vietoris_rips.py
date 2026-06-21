import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
import gudhi as ghd 
from gudhi import plot_persistence_diagram
import matplotlib.ticker as ticker
import random 


os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")

neuron_img = Image.open("mouse_neuron_preprocessed_for_SP.jpg").convert("L")
pixels = np.array(neuron_img)
binary = pixels > 0 

N = 2000 
dimx = binary.shape[0]
dimy = binary.shape[1]
sentinel = np.zeros((dimx, dimy))
pts = []
count = 0 

while count < N:
    pt = [random.randint(0, dimx-1), random.randint(0, dimy-1)]
    if binary[pt[0], pt[1]] == 1:
        pts.append(pt)
        sentinel[pt[0], pt[1]] = 1
        count += 1

plt.imshow(sentinel, cmap="gray")
plt.show()

rips_complex = ghd.RipsComplex(points = pts, max_edge_length=70.0) 
simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)

simplex_tree.compute_persistence()
output = simplex_tree.persistence()
ax = plot_persistence_diagram(output)
ax.set_title("Vietoris-Rips Mouse Neuron", fontsize = 20)
ax.set_aspect("equal")
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel("Birth", fontsize = 18)
plt.ylabel("Death", fontsize = 18)
plt.show() 

#betti curve code
def betti_curve(pairs, dim, r_values):
    relevant = [(birth, death) for d, (birth, death) in pairs if d == dim]
    return [sum(1 for birth, death in relevant 
                if birth <= r and (death > r or death == float('inf'))) for r in r_values]

r_values = [i for i in range(80)]



#plt.plot(r_values, betti_curve(output, 0, r_values), label='H0')
#plt.plot(r_values, betti_curve(output, 1, r_values), label='H1')
#plt.xlabel("Filtration parameter")
#plt.ylabel("Number of generators")
#plt.legend()
#plt.show()