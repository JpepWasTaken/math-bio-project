import numpy as np
import gudhi as ghd
import math
from display import *


# INPUT:  an np.array, a pair of integers (i,j)
# OUTPUT: true if the point with coordinates (i,j) is inside the array
def inside_border(array, i, j):
    return (0 <= i < array.shape[0]) and (0 <= j < array.shape[1])


# INPUT:  an np.array, a triplet of integers (i,j,n)
# OUTPUT: a list of coordinates (x,y) of the points in the n-th Moore neighbourhood around (i,j)
# Ex: the first ring is the 8 pixels surrouding (i,j), etc.
def nth_ring(array, i, j, n):
    res = []
    for k in range(-n, n + 1):
        for l in range(-n, n + 1):
            if (abs(k) == n or abs(l) == n) and inside_border(array, i + k, j + l):
                res.append((i + k, j + l))
    return res


# INPUT:  4 integers corresponding to points (k,l) and (m,n) (in some unspecified array)
# OUTPUT: True if (k,l) and (m,n) are Moore neighbours, i.e. if one is on the 8-square grid around the other.
def are_moore_neighbours(k, l, m, n):
    return max(abs(k - m), abs(l - n)) == 1


# INPUT:  an np.array and a list of 2 integers
# OUTPUT: an integer. Each pixel is mapped to a unique integer, obtained by moving left to right
#         and top to bottom through the array.
# Ex: for a 3 by 3 matrix this returns:
# 0 3 6
# 1 4 7
# 2 5 8
def pixel_to_integer(array, pixel):
    return (pixel[1] * array.shape[0]) + (pixel[0] + 1)


# INPUT:  an array and an integer
# OUTPUT: a pair (i,j)
# This is the inverse function from above. The pair (i,j) is the unique one that gets mapped to n.
def integer_to_pixel(array, n):
    return [n % (array.shape[0]), (n // array.shape[0])]


# INPUT:  a vector [x,y]
# OUTPUT: the square of its norm x^2 + y^2
def square_norm(vec):
    return vec[0] ** 2 + vec[1] ** 2


# INPUT:  a vector [x,y]
# OUTPUT: its principal argument between -pi and pi, in radians
def angle(vec):
    return math.atan2(vec[0], vec[1])


# INPUT:  a pixel [i,j] and a vector [x,y]
# OUTPUT: the filtration value of the pixel in the sweeping plane directed by the vector
# Quirk: the vector is a vector in R^2, so x is abscissa and y is the ordinate.
#        But np.arrays use matrix notation where i is the row number, so it is an ordinate
#        This is why pixel[1] is matched with vec[0] and not vec[1]
# Quirk: the vertical axis is upside down for the np.array (going down increases i) so the
#        2nd term has a -
def get_filtration_value(pixel, vec):
    return pixel[1] * vec[0] - pixel[0] * vec[1]


# INPUT:  a binary array and a vector [x,y]
# OUTPUT: a triplet (res, min_val, max_val)
#         1) res is a matrix of the same shape as array. If array[i,j] == True, the entry of res at row i and column j
#         is the filtration value of the pixel (i,j) in the sweeping plane directed by vec. Otherwise, res[i,j] = infinity.
#         2) min_val is the minimum filtration value of pixels in array in the sweeping plane directed by vec.
#         3) max_val is the maximum filtration value of pixels in array in the sweeping plane directed by vec.
def get_all_filtrations(array, vec):
    ymax, xmax = array.shape[0], array.shape[1]
    min_val, max_val = math.inf, -math.inf
    res = np.zeros((ymax, xmax))
    for j in range(ymax):
        for i in range(xmax):
            if array[j, i] != 0:
                res[j, i] = get_filtration_value([j, i], vec)
                if res[j, i] < min_val:
                    min_val = res[j, i]
                elif res[j, i] > max_val:
                    max_val = res[j, i]
            else:
                res[j, i] = math.inf
    return (res, min_val, max_val)


# INPUT:  a gudhi simplex tree object st, assumed to be empty, a binary array, and a vector
# OUTPUT: nothing. The full simplicial tree of the sweeping plane filtration of array directed by vec
#         is stored in st at the end of the function call.
def construct_simplex_tree(st, array, vec):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if array[i, j] == True:
                # insert 0 simplices
                st.insert(
                    [pixel_to_integer(array, [i, j])],
                    filtration=get_filtration_value([i, j], vec),
                )
                ring = nth_ring(array, i, j, 1)
                # insert 1 simplices
                for k, l in ring:
                    if array[k, l] == True:
                        st.insert(
                            [
                                pixel_to_integer(array, [i, j]),
                                pixel_to_integer(array, [k, l]),
                            ],
                            filtration=max(
                                get_filtration_value([i, j], vec),
                                get_filtration_value([k, l], vec),
                            ),
                        )
                        # insert 2 simplices
                        for m, n in ring:
                            if array[m, n] == True:
                                if are_moore_neighbours(k, l, m, n):
                                    st.insert(
                                        [
                                            pixel_to_integer(array, [i, j]),
                                            pixel_to_integer(array, [k, l]),
                                            pixel_to_integer(array, [m, n]),
                                        ],
                                        filtration=max(
                                            get_filtration_value([i, j], vec),
                                            get_filtration_value([k, l], vec),
                                            get_filtration_value([m, n], vec),
                                        ),
                                    )
    st.make_filtration_non_decreasing()


# INPUT:  a parameter value (float or integer) r and a list, representing the (birth, death) pairs of a barcode decomposition
# OUTPUT: the number of bars in the decomposition that are alive at parameter value r
def count_living_bars(r, birth_death_list):
    return sum(
        [
            1
            for (birth, death) in birth_death_list
            if ((birth <= r) and ((death > r) or death == float("inf")))
        ]
    )


# INPUT:  a list of integers list_
# OUTPUT: the total variation of that list (an integer)
# Ex: the list [1,3,1] will output 1+2+2 = 5.
# The list is meant to be the output of the count_living_bars function.
def compute_total_variation_of_list(list_):
    res = list_[0]
    for i in range(len(list_) - 1):
        res += abs(list_[i] - list_[i + 1])
    return res


# INPUT:  a binary array and a vector [x,y]
# OUTPUT: the total variation of the function f mapping a filtration parameter r to the number of connected components alive at r
#         (in the sweeping plane on array directed by vec)
def total_variation_of_connected_components(array, vec):
    st = ghd.SimplexTree()
    construct_simplex_tree(st, array, vec)
    st.compute_persistence()
    filtration_values = [filtration for simplex, filtration in st.get_filtration()]
    filtration_values.sort()
    # isolate the birth, pair pairs of 0-dimensional simplices
    zeroth_homology_barcode = [
        (birth, death) for dim, (birth, death) in st.persistence() if dim == 0
    ]
    counts = [count_living_bars(r, zeroth_homology_barcode) for r in filtration_values]
    res = compute_total_variation_of_list(counts)
    st.clear()
    return res


# INPUT:  an integer
# OUTPUT: a list of lists [[cos(x), sin(x), x] , [cos(y), sin(y), y]]...
#         The pair (cos(x), sin(x)) is the coordinates of a vector on the unit circle, with angle (in degrees) x
# The vectors are chosen so that (1, 0) and (-1, 0) are always in the list, and the rest are evenly spaced in-between.
def get_vectors(num_vectors):
    res = []
    step = 180 / (num_vectors - 1)
    for i in range(num_vectors):
        angle_deg = i * step
        angle_rad = np.deg2rad(angle_deg)
        res.append([math.cos(angle_rad), math.sin(angle_rad), angle_deg])
    return res


# BUGGY!
# INPUT:  a binary array and an integer
# OUTPUT: a list of lists [[total variation of connected components for direction v, angle of v in degrees], [...]]...
# The total variation (of connected components) is computed for the sweeping planes for each of the vectors
# produced by the get_vectors function. We keep track of the angle of each vector for the wedge display.
def compute_multivector_TV_for_wedge_display(array, num_vectors):
    res = []
    vectors_list = get_vectors(num_vectors)
    count = len(vectors_list)
    curr = 1
    for vec in vectors_list:
        print("Dealing with vector", curr, "out of", count)
        res.append([total_variation_of_connected_components(array, vec), vec[2]])
        curr += 1
    return res


######## Somewhat deprecated
# We can definitely replace this by a function that computes multivector TV for an arbitrary array of vectors


# INPUT:  a binary array and a(n) (empty) simplex tree
# OUTPUT: total variations of connected components in each cardinal and diagonal direction.
def cardinal_and_diagonal_total_variations(st, array):
    res = []
    # iterate through the various sweeping plane directions, computing the total variation for each
    for method_index in range(len(methods_vect)):
        # the max filtration value depends on the specific plane direction chosen
        construct_simplex_tree(st, array, methods_vect[method_index])
        st.compute_persistence()
        h0_gens = [
            (birth, death) for dim, (birth, death) in st.persistence() if dim == 0
        ]
        if method_index <= 3:
            r_max = array.shape[0]
        else:
            r_max = array.shape[0] + array.shape[1]
        r_values = [i for i in range(r_max + 10)]
        counts = [count_living_bars(r, h0_gens) for r in r_values]
        var = compute_total_variation_of_list(counts)
        print(
            "For method ",
            methods_vect[method_index].__name__,
            "the total variation is",
            var,
        )
        res.append(var)
        st.clear()


######### DEPRECATED
# Since we can just compute filtration values omnidirectionally, these specific functions are no longer needed.


# Filtration values of a pixel in 8 directions
def top(array, pixel):
    return pixel[1]


def bot(array, pixel):
    return array.shape[1] - pixel[1]


def left(array, pixel):
    return pixel[0]


def right(array, pixel):
    return array.shape[0] - pixel[0]


def tl_br(array, pixel):
    return pixel[0] + pixel[1]


def tr_bl(array, pixel):
    return abs(array.shape[0] - 1 - pixel[0]) + pixel[1]


def bl_tr(array, pixel):
    return pixel[0] + abs(array.shape[1] - 1 - pixel[1])


def br_tl(array, pixel):
    return abs(array.shape[0] - 1 - pixel[0]) + abs(array.shape[1] - 1 - pixel[1])


methods_vect = [top, bot, left, right, tl_br, tr_bl, bl_tr, br_tl]


def betti_curve(pairs, dim, r_values):
    relevant = [(birth, death) for d, (birth, death) in pairs if d == dim]
    return [
        sum(
            1
            for birth, death in relevant
            if birth <= r and (death > r or death == float("inf"))
        )
        for r in r_values
    ]


### Old stuff for computing sweeping planes on the neurons.
# Commented out in case the examples are still useful.


# Sweeping plane filtration computation for my lad AUSTIN
# Image = Image.open("mouse_neuron_preprocessed_for_SP.jpg").convert("L")
# binary mask!!!!!!! key
# neuron_img = np.array(Image) > 40

# simptree = ghd.SimplexTree()
# construct_simplex_tree(simptree, neuron_img, tl_br)
# simptree.compute_persistence()
# output = simptree.persistence()
# ax = plot_persistence_diagram(output)
# ax.set_title("Top-left/bot-right Sweeping Plane for Neuron", fontsize = 20)
# ax.set_aspect("equal")
# plt.tick_params(axis='both', which='major', labelsize=16)
# plt.xlabel("Birth", fontsize = 18)
# plt.ylabel("Death", fontsize = 18)
# plt.show()

# simptree2 = ghd.SimplexTree()
# construct_simplex_tree(simptree2, neuron_img, left)
# simptree2.compute_persistence()
# output2 = simptree2.persistence()


# r_values = [i for i in range(neuron_img.shape[0] + neuron_img.shape[1])]


# plt.plot(r_values, betti_curve(output, 0, r_values), label='H0')
# plt.plot(r_values, betti_curve(output, 1, r_values), label='H1')
# plt.xlabel("Filtration parameter")
# plt.ylabel("Number of generators")
# plt.legend()
# plt.show()


############# EXAMPLE USE #######################################


# pixels = np.zeros((5,5))
# ring = nth_ring(pixels, 2, 2, 1)
# for k,l in ring:
#    pixels[k,l] = 1
# This is for a 5 by 5 grid with a ring in the middle. Not very interesting - each method has total variation 1.
# But this is the syntax to use - just plug in your own array instead of pixels.


# st = ghd.SimplexTree()
# cardinal_and_diagonal_total_variations(st, pixels)

# Code that kills small connected components of the skeleton image
# labeled = label(pixels)
# for region in regionprops(labeled):
#    if region.area < 40:
#        for i,j in region.coords:
#            pixels[i,j] = False
