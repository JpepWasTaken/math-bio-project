import os 
import skimage as ski 
from PIL import Image 
from skimage.measure import label, regionprops
import numpy as np 
import matplotlib.pyplot as plt
from gudhi import CubicalComplex 
from gudhi import plot_persistence_barcode
from PIL import Image, ImageDraw, ImageFont

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images")


neuron_img = Image.open("mouse_neuron_1_alt.jpg").convert("L")
pixels = np.array(neuron_img) > 0 

b_max = pixels.max() 


def cull_below_threshold(img, thr): 
    dimx, dimy = img.shape
    for i in range(dimx): 
        for j in range(dimy): 
            if img[i,j] < thr: 
                img[i,j] = 0



#prototype of what one iteration of the gif-making process will look like 

#cull_below_threshold(pixels, 50) #kill pixels with brightness lower than 50 
#img = Image.fromarray((pixels * 255).astype(np.uint8), mode='L')

#turn it into a python image, draw text on it, then convert back to an array for display 

#plt.imshow(new_pixels, cmap = "gray")
#plt.show()





#make an animation of the neuron filtrations 
font = ImageFont.truetype("arial.ttf", size=30)
img = Image.fromarray((pixels * 255).astype(np.uint8), mode='L') 
frames = []

#make the filtration parameter go from 0 all the way to 159. Each value, save an image of it. 
for b in range(0,159): 
    neuron_img = np.array(Image.open("mouse_neuron_1_alt.jpg").convert("L"))
    neuron_img = neuron_img >= b
    frame = Image.fromarray((neuron_img * 255).astype(np.uint8), mode='L')
    draw = ImageDraw.Draw(frame)
    draw.text((440, 70), f"r = {159 - b}", fill=255, font=font)  # fill=255 for white text
    frames.append(frame)  

frames.reverse()

os.chdir(r"C:\Users\eudes\Documents\Purdue\Math Bio (Spring 2026)\HW3\images\animation")

for i in range(len(frames)): 
    frames[i].save(f"frame{i}.png")



#frames[0].save(
#    'neuron_growth.png',
#    save_all=True,
#    append_images=frames[1:],
#    duration=200,    # milliseconds per frame
#    loop=0          # 0 = loop forever
#)