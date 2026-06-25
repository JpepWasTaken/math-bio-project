import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
import gudhi as ghd 
from plane_filtrations import * 
from display import * 

os.chdir("images")




pixels = np.zeros((5,5))
ring = nth_ring(pixels, 2, 2, 1)
for k,l in ring:
    pixels[k,l] = 1

retina_image = Image.open("nardini_retina.png").convert("L")
binarized_retina_image = ski.morphology.skeletonize(np.array(retina_image) > 100) 



#create_sweeping_plane_gif(binarized_retina_image, [1,0], "left_to_right.gif")


#plt.imshow(pixels, cmap='gray')
#plt.show()
#plt.clf() 

results = compute_tv_for_display(binarized_retina_image, 100)
compute_wedge_diagram(results)



#simptree = ghd.SimplexTree() 
#insert_simplices(simptree, pixels, [1,2])
#simptree.compute_persistence()
#output = simptree.persistence()
#ax = plot_persistence_diagram(output)
#ax.set_title("Left to Right sweeping plane", fontsize = 20)
#ax.set_aspect("equal")
#plt.tick_params(axis='both', which='major', labelsize=16)
#plt.xlabel("Birth", fontsize = 18)
#plt.ylabel("Death", fontsize = 18)
#plt.show() 

vals = np.zeros((5,5))
for i in range(5):
    for j in range(5): 
        vals[i,j] = pixel_to_integer(pixels, [i,j])



#st = ghd.SimplexTree()
#mp.insert_simplices(st, pixels, mp.top)
#print(st.num_simplices())

#for l,f in st.get_simplices(): 
#    print("Simplex found:", l, "with filtration value", f)


#st.compute_persistence()
#pers_list = st.persistence() 

#ghd.plot_persistence_diagram(st.persistence())
#plt.show()