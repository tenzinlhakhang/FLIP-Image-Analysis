%matplotlib inline
import os
from matplotlib import pyplot as plt
import mahotas
import numpy as np
from scipy import ndimage
from skimage import feature
from skimage.draw import circle_perimeter
from skimage.util import img_as_uint, img_as_float
from IPython.html import widgets
from skimage import filter
from skimage.filter import threshold_adaptive
from skimage.morphology import (closing, opening, selem, remove_small_objects, erosion, label, watershed, binary_dilation, black_tophat)
from skimage.measure import perimeter, regionprops
from skimage.segmentation import clear_border
from skimage.exposure import (equalize_hist, histogram)
from skimage import morphology
from skimage import img_as_bool
from skimage.filter.rank import median, gradient, enhance_contrast
from skimage.morphology import disk
from skimage.filter.rank import median
from skimage import restoration
from operator import truediv
from skimage.morphology import skeletonize
from skimage.draw import circle



def odd_number(x):
    if x % 2 == 1:
        return True
    else:
        return False


# Functions for naturally sorting filenames

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]

# Function to return coordinates of unusually bright debris/beads that may
# be skewing the quantification of the fluorescence.

def bright_check_ab(image,x,y):
    original_image = mahotas.imread(image)
    filtered = restoration.denoise_tv_bregman(original_image,10)
    struct_element = selem.diamond(3)
    smooth = filter.gaussian_filter(filtered,sigma=0.5)
    adapted = threshold_adaptive(smooth,block_size=300)
    eroded = erosion(adapted,struct_element)
    boolean = img_as_bool(eroded)
    boolsum = np.sum(boolean)
    filled = ndimage.binary_fill_holes(boolean)
    distance = ndimage.distance_transform_edt(filled)
    local_maxi = feature.peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                                labels=filled)
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(distance, markers, mask=filled)
    props = regionprops(labels,original_image)

    coordinates_remove = []
    max_intensity = []
    mean_intensity = []
    min_intensity = []

    for labels in props:
        max_intensity.append(labels.max_intensity)
        mean_intensity.append(labels.mean_intensity)
        min_intensity.append(labels.min_intensity)

        if labels.max_intensity > x:
            coordinates_remove.append(labels.coords)
        if labels.max_intensity > y:
            coordinates_remove.append(labels.coords)
            
    coord_ = []
    coord_2 = []

    for inner_l in coordinates_remove:
        for item in inner_l:
            coord_.append(item)

    for inner_l in coord_:
        for item in inner_l:
            coord_2.append(item)

    y_coordinates = z[1::2]
    x_coordinates = z[0::2]
    
    if len(coordinates_remove) < 8:
        return x_coordinates, y_coordinates
    else:
        return [1,1],[1,1]

# Creates multiple plots of the original image displaying histogram intensity values for
# the regions of interest.

def plot(original,ab_name,ring_bead_intensity, background_intensity, bead_image, background_image,bins1,bins2):
    
    fig,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(nrows=1,ncols=5,figsize=(18,5))
    
    ax1.imshow(original,cmap=plt.cm.gray, interpolation='none') 
    ax1.set_title(ab_name,fontsize=24)
    ax1.axis('off')

    ax2.hist(ring_bead_intensity,bins=bins1)
    ax2.set_title('Control_Outter_Ring',fontsize=24)

    ax3.hist(background_intensity,bins=bins2)
    ax3.set_title('Control_Inner_Bead',fontsize=24)

    ax4.imshow(bead_image,cmap=plt.cm.gray,interpolation='none')
    ax4.set_title('Bead_signal',fontsize=24)
    ax4.axis('off')

    ax5.imshow(background_image,cmap=plt.cm.gray,interpolation='none')
    ax5.set_title('Background_signal',fontsize=24)
    ax5.axis('off')
    return fig
    
# returns unzipped x,y coordinates from centroid coordinates[0,1].
def xycoord(coord):
        
    bead_coord = []
    bead_coord_two = []
    
    for inner_l in coord:
        for item in inner_l:
            bead_coord.append(item)

    for inner_l in bead_coord:
        for item in inner_l:
            bead_coord_two.append(item)

    y_coordinates = bead_coord_two[1::2]
    x_coordinates = bead_coord_two[0::2]   
    
    return x_coordinates, y_coordinates

# Labels the input image using the intensity values from the original image.
# returns the properties of the labeled objects in the input image.


def labeling(image,original_image):

    distance = ndimage.distance_transform_edt(image)
    local_maxi = feature.peak_local_max(distance, indices=False, footprint=np.ones((1, 1)),
                                labels=image)
    markers = ndimage.label(local_maxi)[0]
    labels = watershed(distance, markers, mask=image)
    props = regionprops(labels,original_image)
    return props

# Used to remove x,y coordinates in the input image.

def remove_coordinates(image,x_coordinates,y_coordinates):
    for val in x_coordinates, y_coordinates:
        image[x_coordinates,y_coordinates] = False
        
# Function for determining fluorescence value for the control image
def control_function(image):
    #coordinates to remove
    x_coordinates,y_coordinates = bright_check_control(image)  
    name_string = str(image)
    original_image = mahotas.imread(image)
    denoised = restoration.denoise_tv_chambolle(original_image,9)
    adapt_threshold = threshold_adaptive(denoised,block_size=150)
    struct_element = selem.diamond(2)
    opened = morphology.opening(adapt_threshold,struct_element)
    bead_value = ndimage.binary_fill_holes(opened)
    removed = remove_small_objects(bead_value,min_size=400)

    beads = clear_border(removed)

    # Removes x,y coordinates which were found to be debris/anomaly.
    remove_coordinates(beads,x_coordinates,y_coordinates)
    
	# Labels the beads in the image without the bright debris.
    beads_label = labeling(beads,original_image)
    centroid_to_make_circle = []
    diameter_beads_remaining = []
    beads_to_remove = [] 
	
	# If beads labeled pass an area/circularity threshold then
	# the centroid and diameter of the bead is appended to a list.
	# If threshold not passed, the coordinates of those beads are
	# appended to a list to be removed.
	
    for labels in beads_label:
        if labels.area > 900 and labels.eccentricity <0.7:
            centroid_to_make_circle.append(labels.centroid)
            diameter_beads_remaining.append(labels.equivalent_diameter)
        else:
            beads_to_remove.append(labels.coords)
    x_coordinates, y_coordinates = xycoord(beads_to_remove)
    remove_coordinates(beads,x_coordinates,y_coordinates)

	# Unzip the centroid coordinate for later use
    centroid_one = []
    for inner_l in centroid_to_make_circle:
        for item in inner_l:
            centroid_one.append(item)
    y_circle_coordinates = centroid_one[1::2]
    x_circle_coordinates = centroid_one[0::2]       

    # Determine the shape of the original image
    # So as it create numpy array of same shape for later use.
    x,y = original_image.shape
    img = np.zeros((x, y), dtype=np.bool)

	# Find ~radius of each object from the list with object diameter values.
    radius = []
    for diameter in diameter_beads_remaining:
        radius.append(diameter/2.5)
	# Create circles in new img array with centroid coordinates from labeled objects
	# As well as previously calculated radius
	
    for x,y,z in zip(x_circle_coordinates,y_circle_coordinates,radius):
        rr, cc = circle(x,y,z)
        img[rr,cc] = 1

	# Find properties of labeled object using original_image as intensity reference.
    drawn_circle_props = labeling(img,original_image)
    circle_coordinates = []
    mean_intensity_drawncircle = []
    
	# Find coordinates for drawn circles and mean intensity
    for label in drawn_circle_props:
        circle_coordinates.append(label.coords)
        mean_intensity_drawncircle.append(labels.mean_intensity)
        
	# Average intensity values for the drawn in circle inside the bead
    inner_bead_avg = np.mean(mean_intensity_drawncircle)
    # return x,y coordinates from circle [0,1] coordinate.
    x_coordinates, y_coordinates = xycoord(circle_coordinates)
	# Remove x,y coordinates from original beads image.
	# This leaves donut shaped objects. 
	# Isolating the region of interest which is the ring of the bead.
	
    remove_coordinates(beads,x_coordinates,y_coordinates)
    ring_intensity = original_image[beads]
    inner_bead_intensity = original_image[img]
    
    minint = np.min(ring_intensity)
    maxint = np.max(ring_intensity)
    bins1 = np.linspace(minint,maxint,256)

    minint = np.min(inner_bead_intensity)
    maxint = np.max(inner_bead_intensity)
    bins2 = np.linspace(minint,maxint,256)
    
    plot(original_image,name_string, ring_intensity, inner_bead_intensity, beads, img,bins1=bins1,bins2=bins2)
    
    
    igg_control.append(np.mean(ring_intensity)- np.mean(inner_bead_intensity))
    return np.mean(ring_intensity), np.mean(inner_bead_intensity), np.mean(ring_intensity)- np.mean(inner_bead_intensity)
    

def ab_function(image,x,y):

    x_coordinates,y_coordinates = bright_check_ab(image)
    name_string = str(image)
    original_image = mahotas.imread(image)
    adapt_threshold = threshold_adaptive(original_image,block_size=x)
    struct_element = selem.diamond(y)
    opened = morphology.opening(adapt_threshold,struct_element)
    bead_value = ndimage.binary_fill_holes(opened)
    removed = remove_small_objects(bead_value,min_size=500)
    boolean = img_as_bool(removed)
    beads = clear_border(removed)
    remove_coordinates(beads,x_coordinates,y_coordinates)
    #beads to remove
    beads_label = labeling(beads,original_image)
    centroid_to_make_circle = []
    diameter_beads_remaining = []
    beads_to_remove = [] 

    for labels in beads_label:
        if labels.area > 900 and labels.eccentricity <0.7:
            centroid_to_make_circle.append(labels.centroid)
            diameter_beads_remaining.append(labels.equivalent_diameter)
        else:
            beads_to_remove.append(labels.coords)

    x_coordinates, y_coordinates = xycoord(beads_to_remove)
    remove_coordinates(beads,x_coordinates,y_coordinates)

    plt.imshow(beads)
    plt.gray()

    centroid_one = []

    for inner_l in centroid_to_make_circle:
        for item in inner_l:
            centroid_one.append(item)

    y_circle_coordinates = centroid_one[1::2]
    x_circle_coordinates = centroid_one[0::2]       

    x,y = original_image.shape
    img = np.zeros((x, y), dtype=np.bool)

    radius = []

    for diameter in diameter_beads_remaining:
        radius.append(diameter/2.5)

    for x,y,z in zip(x_circle_coordinates,y_circle_coordinates,radius):
        rr, cc = circle(x,y,z)
        img[rr,cc] = 1

    drawn_circle_props = labeling(img,original_image)

    # Finding coordinates for drawn in circles
    circle_coordinates = [] 
    mean_intensity_drawncircle = []
    
    for label in drawn_circle_props:
        circle_coordinates.append(label.coords)
        mean_intensity_drawncircle.append(labels.mean_intensity)
        
    inner_bead_avg = np.mean(mean_intensity_drawncircle)
    x_coordinates, y_coordinates = xycoord(circle_coordinates)
    remove_coordinates(beads,x_coordinates,y_coordinates)

    ring_intensity = original_image[beads]
    inner_bead_intensity = original_image[img]
    
    minint = np.min(ring_intensity)
    maxint = np.max(ring_intensity)
    bins1 = np.linspace(minint,maxint,256)

    minint = np.min(inner_bead_intensity)
    maxint = np.max(inner_bead_intensity)
    bins2 = np.linspace(minint,maxint,256)
    
    plot(original_image,name_string, ring_intensity, inner_bead_intensity, beads, img,bins1=bins1,bins2=bins2)

    igg_control.append(np.mean(ring_intensity)- np.mean(inner_bead_intensity))
    return np.mean(ring_intensity), np.mean(inner_bead_intensity), np.mean(ring_intensity)- np.mean(inner_bead_intensity)


# Function used to isolate and quantify fluorescence in antibody images
# Which passes first filter for determining presence of protein binding.

def FLIP_greater(image):

    name_string = str(image)
    x_coordinates, y_coordinates = bright_check_ab(image)
    original_image = mahotas.imread(image)
    filtered = restoration.denoise_tv_bregman(original_image,.1)
    struct_element = selem.diamond(2)
    smooth = filter.gaussian_filter(filtered,sigma=0.5)
    adapted = threshold_adaptive(smooth,block_size=55)
    adaptedbg = threshold_adaptive(smooth,block_size=60)
    eroded = erosion(adapted,struct_element)
    opened = opening(adapted,struct_element)
    openedbg = opening(adaptedbg,struct_element)
    ndstruct = selem.diamond(1)   
    ndimg = ndimage.binary_opening(eroded,ndstruct)
    skeleton = skeletonize(ndimg)

    plt.imshow(skeleton)
    if np.sum(skeleton) < 1000:
        FLIP_no_skeleton(image)
    else:
        adaptedbg = threshold_adaptive(smooth,block_size=60)
        openedbg = opening(adaptedbg,struct_element)
        dilated = binary_dilation(openedbg,struct_element)
        closing = ndimage.binary_closing(dilated)
        w = ~np.array(closing)
        summed = ndimg

        for val in x_coordinates,y_coordinates:
            w[x_coordinates,y_coordinates] = False

        for val in x_coordinates,y_coordinates:
            skeleton[x_coordinates,y_coordinates] = False

        ring_intensity = original_image[skeleton]
        background_intensity = original_image[w]

        minint = np.min(ring_intensity)
        maxint = np.max(ring_intensity)
        bins1 = np.linspace(minint,maxint,256)

        minint = np.min(background_intensity)
        maxint = np.max(background_intensity)
        bins2 = np.linspace(minint,maxint,256)

        FLIP_ab_value.append(np.mean(ring_intensity)-np.mean(background_intensity))
        plot(original_image,name_string, ring_intensity, background_intensity, skeleton, w,bins1=bins1,bins2=bins2)
        plt.tight_layout()

        return np.mean(ring_intensity), np.mean(background_intensity), np.mean(ring_intensity)-np.mean(background_intensity)

def FLIP_no_skeleton(image):
    
    name_string = str(image)
    x_coordinates, y_coordinates = bright_check_ab(image)
    original_image = mahotas.imread(image)
    filtered = restoration.denoise_tv_bregman(original_image,1)
    struct_element = selem.diamond(1)
    smooth = filter.gaussian_filter(original_image,sigma=1.3)
    adapted = threshold_adaptive(smooth,block_size=100)
    adaptedbg = threshold_adaptive(smooth,block_size=55)
    eroded = erosion(adapted,struct_element)
    opened = opening(adapted,struct_element)
    openedbg = opening(adaptedbg,struct_element)
    ndstruct = selem.diamond(1)
    ndimg = ndimage.binary_opening(eroded,ndstruct)
    skeleton = skeletonize(ndimg)   
    plt.imshow(skeleton)
    adaptedbg = threshold_adaptive(smooth,block_size=60)
    openedbg = opening(adaptedbg,struct_element)
    dilated = binary_dilation(openedbg,struct_element)
    closing = ndimage.binary_closing(dilated)

    w = ~np.array(closing)
    summed = ndimg

    for val in x_coordinates,y_coordinates:
        ndimg[x_coordinates,y_coordinates] = False

    for val in x_coordinates,y_coordinates:
        w[x_coordinates,y_coordinates] = False

    ring_intensity = original_image[ndimg]
    background_intensity = original_image[w]

    minint = np.min(ring_intensity)
    maxint = np.max(ring_intensity)
    bins1 = np.linspace(minint,maxint,256)

    minint = np.min(background_intensity)
    maxint = np.max(background_intensity)
    bins2 = np.linspace(minint,maxint,256)

    FLIP_ab_value.append(np.mean(ring_intensity)-np.mean(background_intensity))
    plot(original_image,name_string, ring_intensity, background_intensity, ndimg, w,bins1=bins1,bins2=bins2)
    print np.mean(ring_intensity), np.mean(background_intensity), np.mean(ring_intensity)-np.mean(background_intensity)
    plt.tight_layout()
