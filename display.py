import matplotlib.pyplot as plt 
import numpy as np 
import random as rd
from math import pi 

# idea: take as input an array vec = [[value, angle], [value2, angle2], [value3, angle3]]
# and display it as a half-circle. angles are expected to be in degrees 

def compute_wedge_diagram(vec):
    vec.sort(key = lambda x : x[1])  #sorts the array based on the angle values 
    angles_list = np.array([d[1] for d in vec])
    val_list = np.array([d[0] for d in vec])
    step = angles_list[1] - angles_list[0]
    wedge_bound_angles = np.deg2rad(np.array([angles_list[0] - step/2 + i * step for i in range(len(angles_list)+1)]))

    print("angles:", angles_list)
    print("values:", val_list)
    print("bounds:", wedge_bound_angles)
    r_inner, r_outer = 9, 10
    radius_edges = [r_inner, r_outer]

    # this mesh construction expects values in radians 
    Theta, R = np.meshgrid(wedge_bound_angles, radius_edges)
    print("THETA", Theta)
    print(R)
    values_2d = val_list.reshape(1, -1)
    
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(projection='polar')
    mesh = ax.pcolormesh(Theta, R, values_2d, cmap='plasma', shading='flat',
                      edgecolors='white', linewidth=0.5)
    ax.set_thetamin(angles_list[0])
    ax.set_thetamax(angles_list[-1])
    ax.set_ylim(r_inner-2, r_outer)
    fig.colorbar(mesh, ax=ax, pad=0.1, shrink=0.7)
    plt.tight_layout()
    plt.grid(False)
    plt.show()


#step = 180 / 100
#vec = [[rd.uniform(0.5, 1.5), i * step] for i in range(101)]
#compute_wedge_diagram(vec)