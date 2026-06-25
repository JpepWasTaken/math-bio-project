import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
import gudhi as ghd 
from gudhi import plot_persistence_diagram
import math 
from display import *

#checks if an index i,j lies within an array or outside of it 
def out_of_border(array, i, j):
    xmax, ymax = array.shape 
    if i < 0 or i >= xmax or j < 0 or j >= ymax:
        return True
    return False

#takes as an input the coordinates of a pixel (i,j) and an integer n, and returns the coordinate of the 
#pixels in the nth ring around (i,j). The first ring is the 8 pixels surrouding (i,j), etc. 
def nth_ring(array,i,j,n): 
    res = []
    for k in range(-n, n+1):
        for l in range(-n, n+1):
            if (abs(k) == n or  abs(l) == n) and not out_of_border(array, i+k, j+l): 
                #we are on the border of the nth ring, and not oob
                res.append((i+k, j+l))
    return res

#Code that maps pixels on a grid to an integer one-to-one. Goes column by column, starting at the top.  
def pixel_to_integer(array, pixel):
    return pixel[1]*array.shape[0] + pixel[0]

#Recovers coordinates from the integers. 
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

# General helper functions 
def square_norm(vec):
    return (vec[0]**2 + vec[1]**2)

def angle(vec):
    return math.atan2(vec[0], vec[1])

# gets filtration value for an arbitrary pixel + vector 
# the filtration value is just the dot product with the vector, but need to be careful because the 
# y-axis is upside-down in arrays
def get_filtration_value(pixel, vec): 
    py, px = pixel[0], pixel[1] 
    return (px * vec[0] - py * vec[1])

def get_all_filtrations(array, vec): 
    ymax, xmax = array.shape[0], array.shape[1]
    min_val = math.inf 
    res = np.zeros((ymax, xmax))
    for j in range(ymax):
        for i in range(xmax):
            if array[j,i] != 0: 
                res[j,i] = get_filtration_value([j, i], vec)
                if res[j,i] < min_val:
                    # track the lowest filtration value to make all filtration values positive at the end 
                    min_val = res[j,i]
            else: 
                res[j,i] = math.inf 
    print(min_val)
    res = res - min_val * np.ones((ymax, xmax))
    return res 


methods_vect = [top, bot, left, right, tl_br, tr_bl, bl_tr, br_tl]

#Checks if 2 pixels are Moore neighbours (i.e. if one is on the 8-pixel ring around the other)
def are_moore_neighbours(k, l, m, n):
    return (max(abs(k-m), abs(l-n)) == 1)

#Builds the simplex tree object from a binary array then adds simplices to it, using filtration values obtained by 
#applying the "method" function. method is intended to be one of the 8 functions above, corresponding to sweeping 
#planes from the top, bottom, right, left, and the 4 diagonals. 
def insert_simplices(st, array, vec):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if array[i,j] == True: 
                #insert 0 simplices first
                st.insert([pixel_to_integer(array, [i,j])], filtration = get_filtration_value([i,j], vec))
                ring = nth_ring(array, i, j, 1)
                #insert 2 simplices next
                for (k,l) in ring: 
                    if array[k,l] == True: 
                        st.insert([pixel_to_integer(array, [i,j]), pixel_to_integer(array, [k,l])], 
                                  filtration = max(get_filtration_value([i,j], vec), get_filtration_value([k,l], vec)))
                        for (m,n) in ring: 
                            if array[m,n] == True:
                                if are_moore_neighbours(k,l,m,n): 
                                    st.insert([pixel_to_integer(array, [i,j]), 
                                               pixel_to_integer(array, [k,l]), pixel_to_integer(array, [m,n])], 
                                               filtration = max(get_filtration_value([i,j], vec), get_filtration_value([k,l], vec), 
                                                                get_filtration_value([m,n], vec)))
    st.make_filtration_non_decreasing()


#Helper function. For a given parameter value r and h0_gens the list of dim 0 (birth, death) pairs, computes how many dim 0 
#bars are alive at r
def count_h0_gens(r, h0_gens):
    return sum([1 for (birth,death) in h0_gens if ((birth <= r) and ((death > r) or death == float("inf")))])   


#Helper function. Computes total variation of a vector. Intended to take the output of the above function as its 
#input. 
def total_var(h0_counts):
    res = h0_counts[0]
    for i in range(len(h0_counts)-1):
        res += abs(h0_counts[i] - h0_counts[i+1])
    return res


def single_plane_0_dim_total_variation(array, vec):
    st = ghd.SimplexTree()
    insert_simplices(st, array, vec)
    st.compute_persistence()
    filtration_values = [filtration for simplex, filtration in st.get_filtration()]
    filtration_values.sort() 
    # isolate the birth, pair pairs of 0-dimensional simplices 
    h0_gens = [(birth,death) for dim, (birth,death) in st.persistence() if dim == 0] 
    counts = [count_h0_gens(r, h0_gens) for r in filtration_values]
    var = total_var(counts)
    st.clear() 
    return var 


def get_vectors(num_vectors):
    res = []
    step = 180/(num_vectors-1)
    for i in range(num_vectors):
        angle_deg = i * step 
        angle_rad = np.deg2rad(angle_deg)
        res.append([math.cos(angle_rad), math.sin(angle_rad), angle_deg])
    return res


def compute_tv_for_display(array, num_vectors):
    res = [] # array of [TV in direction v, angle of v]
    vectors_list = get_vectors(num_vectors)
    for vec in vectors_list: 
        res.append([single_plane_0_dim_total_variation(array, vec), vec[2]])
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
        print("For method ", methods_vect[method_index].__name__, "the total variation is", var)
        res.append(var)
        st.clear()


def betti_curve(pairs, dim, r_values):
    relevant = [(birth, death) for d, (birth, death) in pairs if d == dim]
    return [sum(1 for birth, death in relevant 
                if birth <= r and (death > r or death == float('inf'))) for r in r_values]










### Old stuff for computing sweeping planes on the neurons. 
# Commented out in case the examples are still useful. 



#Sweeping plane filtration computation for my lad AUSTIN 
#Image = Image.open("mouse_neuron_preprocessed_for_SP.jpg").convert("L")
#binary mask!!!!!!! key 
#neuron_img = np.array(Image) > 40

#simptree = ghd.SimplexTree() 
#insert_simplices(simptree, neuron_img, tl_br)
#simptree.compute_persistence()
#output = simptree.persistence()
#ax = plot_persistence_diagram(output)
#ax.set_title("Top-left/bot-right Sweeping Plane for Neuron", fontsize = 20)
#ax.set_aspect("equal")
#plt.tick_params(axis='both', which='major', labelsize=16)
#plt.xlabel("Birth", fontsize = 18)
#plt.ylabel("Death", fontsize = 18)
#plt.show() 

#simptree2 = ghd.SimplexTree() 
#insert_simplices(simptree2, neuron_img, left)
#simptree2.compute_persistence()
#output2 = simptree2.persistence()


#r_values = [i for i in range(neuron_img.shape[0] + neuron_img.shape[1])]



#plt.plot(r_values, betti_curve(output, 0, r_values), label='H0')
#plt.plot(r_values, betti_curve(output, 1, r_values), label='H1')
#plt.xlabel("Filtration parameter")
#plt.ylabel("Number of generators")
#plt.legend()
#plt.show()












############# EXAMPLE USE #######################################


#pixels = np.zeros((5,5))
#ring = nth_ring(pixels, 2, 2, 1)
#for k,l in ring:
#    pixels[k,l] = 1
#This is for a 5 by 5 grid with a ring in the middle. Not very interesting - each method has total variation 1. 
#But this is the syntax to use - just plug in your own array instead of pixels. 


#st = ghd.SimplexTree()
#total_variations(st, pixels)

#Code that kills small connected components of the skeleton image
#labeled = label(pixels)
#for region in regionprops(labeled):
#    if region.area < 40: 
#        for i,j in region.coords:
#            pixels[i,j] = False