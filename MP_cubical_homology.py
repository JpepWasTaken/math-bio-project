import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
from gudhi import CubicalComplex 
from gudhi import plot_persistence_diagram

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")

#Computes the cubical persistence associated with the neuron image. 
#TODO: add images of what the neuron looks like with a variety of brightness parameters. Will help to interpret data. 



neuron_img = Image.open("mouse_neuron_1_alt.jpg").convert("L")

pixels = np.array(neuron_img)


cc = CubicalComplex(top_dimensional_cells= pixels.max() - pixels)
output = cc.persistence()

#plot_persistence_diagram(cc.persistence(), max_intervals=1000000)
#plt.show()
#plt.xlabel("Parameter r")


def betti_curve(pairs, dim, r_values):
    relevant = [(birth, death) for d, (birth, death) in pairs if d == dim]
    return [sum(1 for birth, death in relevant 
                if birth <= r and (death > r or death == float('inf'))) for r in r_values]

r_values = [i for i in range(170)]


plt.plot(r_values, betti_curve(output, 0, r_values), label='H0')
plt.plot(r_values, betti_curve(output, 1, r_values), label='H1')
plt.xlabel("Filtration parameter")
plt.ylabel("Number of generators")
plt.legend()
plt.show()




#txt = "The parameter r is inverse brightness: bright pixels enter the filtration early and dim ones late. Intervals shorter than 3 are culled."
#plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=12)
#plt.show()

