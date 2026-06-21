import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
import gudhi as ghd 
import cv2 

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")


preprocessed_image = Image.open("wing1.png").convert("L")
#binary mask!!!!!!! key 
prepped = np.array(preprocessed_image) > 40


colored_neuron_img = Image.open("mouse_neuron_1_alt.jpg")
neuron_img = Image.open("mouse_neuron_1_alt.jpg").convert("L")

colored_pixels = np.array(colored_neuron_img)
pixels = np.array(neuron_img)

thresh = cv2.threshold(pixels, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[0]
print(thresh)
binary = pixels > 40 
#plt.imshow(binary, cmap="gray")
#plt.show()


skeleton = ski.morphology.skeletonize(prepped)




def out_of_border(array, i, j):
    xmax, ymax = array.shape 
    if i < 0 or i >= xmax or j < 0 or j >= ymax:
        return True
    return False

#this function takes as an input the coordinates of a pixel (i,j) and an integer n, and returns the coordinate of the 
#pixels in the nth ring around (i,j). The first ring is the 8 pixels surrouding (i,j), etc. 
def nth_ring(array,i,j,n): 
    res = []
    for k in range(-n, n+1):
        for l in range(-n, n+1):
            if (abs(k) == n or  abs(l) == n) and not out_of_border(array, i+k, j+l): 
                #we are on the border of the nth ring, and not oob
                res.append((i+k, j+l))
    return res


def cull_noise (array):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if array[i,j] == 1: 
                for n in range(1,2):
                    #check the 1st, 2nd, 3rd and 4th rings around i,j. If one of them is entirely 0, cull i,j. 
                    ring = nth_ring(array, i, j, n)
                    if sum([array[k,l] for k,l in ring]) <= 1: 
                        array[i,j] = False 
                        pass 

cull_noise(skeleton)

#Code that kills small connected components of the skeleton image
labeled = label(skeleton)
for region in regionprops(labeled):
    if region.area < 40: 
        for i,j in region.coords:
            skeleton[i,j] = False

#final_no_noise = Image.fromarray(skeleton.astype(np.uint8) * 255)
#final_no_noise.save("mouse_neuron_preprocessed_for_SP.jpg")

#Code that displays overlayed images 
plt.imshow(skeleton, cmap = "gray")
#plt.imshow(colored_pixels, alpha=0.5)
plt.show()
    


#Code that maps pixels on a grid to an integer one-to-one 
def pixel_to_integer(array, pixel):
    return pixel[1]*array.shape[0] + pixel[0]

def integer_to_pixel(array, n): 
    return [n % (array.shape[0]), (n//array.shape[0])]

#Filtration values of a pixel in 8 directions
def top(array, pixel):
    return pixel[1]

def bot(array, pixel):
    return array.shape[1]-pixel[1]
    
def left(array, pixel):
    return pixel[0]

def right(array, pixel):
    return array.shape[0]-pixel[0]

def tl_br(array, pixel): 
    return pixel[0]+pixel[1]

def tr_bl(array, pixel):
    return abs(array.shape[0]-1-pixel[0]) + pixel[1]

def bl_tr(array, pixel):
    return pixel[0] + abs(array.shape[1]-1-pixel[1]) 

def br_tl(array, pixel):
    return abs(array.shape[0]-1-pixel[0]) + abs(array.shape[1]-1-pixel[1]) 

def are_moore_neighbours(k, l, m, n):
    return (max(abs(k-m), abs(l-n)) == 1)

methods_vect = [top, bot, left, right, tl_br, tr_bl, bl_tr, br_tl]

#Code that builds a simplex tree for my sweeping plane 

st = ghd.SimplexTree()

def insert_simplices(st, array, method):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if array[i,j] == True: 
                #insert 0 simplices first
                st.insert([pixel_to_integer(array, [i,j])], filtration = method(array, [i,j]))
                ring = nth_ring(array, i, j, 1)
                #insert 2 simplices next
                for (k,l) in ring: 
                    if array[k,l] == True: 
                        st.insert([pixel_to_integer(array, [i,j]), pixel_to_integer(array, [k,l])], 
                                  filtration = max(method(array, [i,j]), method(array, [k,l])))
                        for (m,n) in ring: 
                            if array[m,n] == True:
                                if are_moore_neighbours(k,l,m,n): 
                                    st.insert([pixel_to_integer(array, [i,j]), 
                                               pixel_to_integer(array, [k,l]), pixel_to_integer(array, [m,n])], 
                                               filtration = max(method(array, [i,j]), method(array, [k,l]), 
                                                                method(array, [m,n])))
    st.make_filtration_non_decreasing()


def count_h0_gens(r, h0_gens):
    return sum([1 for (birth,death) in h0_gens if ((birth <= r) and ((death > r) or death == float("inf")))])    

def total_var(h0_counts):
    res = h0_counts[0]
    for i in range(len(h0_counts)-1):
        res += abs(h0_counts[i] - h0_counts[i+1])
    return res

#Computes the total variation. st is a freshly initialized simplextree object, array the binary picture
def total_variations(st, array):
    res = []
    #iterate through the various sweeping plane directions, computing the total variation for each 
    for method_index in range(len(methods_vect)): 
        #the max filtration value depends on the specific plane direction chosen 
        insert_simplices(st, array, methods_vect[method_index])
        st.compute_persistence()
        h0_gens = [(birth,death) for dim, (birth,death) in st.persistence() if dim == 0] 
        if method_index <= 3:
            r_max = array.shape[0]
        else:
            r_max = array.shape[0] + array.shape[1]
        r_values = [i for i in range(r_max+10)]
        counts = [count_h0_gens(r, h0_gens) for r in r_values]
        var = total_var(counts)
        print("For method", methods_vect[method_index].__name__, "the total variation is", var)
        res.append(var)
        st.clear()
    return res 

res = total_variations(st, skeleton)
print("top-bot", res[0]/res[1])
print("left-right", res[2]/res[3])
print("SE", res[4]/res[7])
print("SW", res[5]/res[6])



### Pretty sure this is useless now 


#insert_simplices(st, skeleton, bl_tr)

#st.compute_persistence()
#pers_list = st.persistence() 
#h0_gens = [(birth,death) for dim, (birth,death) in pers_list if dim == 0] 

#max_r = skeleton.shape[0]+skeleton.shape[1]


#r_values = [i for i in range(max_r+10)]
#counts = [count_h0_gens(r) for r in r_values]

#print(total_var(counts))
#plt.plot(r_values, counts)
#plt.xlabel("filtration parameter")
#plt.ylabel("number of connected components")
#plt.show()