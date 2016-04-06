import os
from flask import Flask, request, redirect, url_for, send_file
from werkzeug import secure_filename
from flask import send_from_directory
from app import app
import matplotlib.pyplot as plt
import numpy as np
from StringIO import StringIO
from flask import render_template
import mahotas
import skimage
from flask import make_response, request
from skimage import data , io
from imread import imread_from_blob
from skimage.filter import threshold_adaptive
import uuid
from scipy import ndimage
from skimage import feature
from skimage.draw import circle_perimeter
from skimage.filter import threshold_adaptive
from skimage.morphology import (closing, opening, selem, remove_small_objects, erosion, label, watershed, binary_dilation, black_tophat)
from skimage.measure import perimeter, regionprops
from skimage.segmentation import clear_border
from skimage import morphology
from skimage import img_as_bool
from skimage.filter.rank import median, gradient, enhance_contrast
from skimage.morphology import disk
from skimage.filter.rank import median
from skimage import restoration
from operator import truediv
from skimage.morphology import skeletonize
from image import *
from flask import g


UPLOAD_FOLDER = '/Users/tclhakhang/Desktop/FLIP Software/tmp'
ALLOWED_EXTENSIONS = set(['tif', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename):
            stored_filename = str(file.filename)
            filename = str(uuid.uuid1()) + '.' + file.filename.rsplit('.',1)[1] + '*' + stored_filename
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return redirect(url_for('images',
                                    filename=filename))
    
    return render_template("homepage.html")
    
@app.route('/images/<filename>')
def images(filename):
	return render_template("analysis.html", title=filename)

 	
	
@app.route('/test/<filename>')
def test(filename):
	image_name = str(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	original_name = filename.split('*')
	original = original_name[1]
	
	figimg = FLIP_greater(image_name)
	
	img = StringIO()
	plt.savefig(img)
	img.seek(0)
	
	return send_file(img, mimetype='image/png')


	
@app.route('/test2/<filename>')
def test2(filename):
	image_name = str(os.path.join(app.config['UPLOAD_FOLDER'], filename))
	original_name = filename.split('*')
	original = original_name[1]
	figimg = FLIP_greater(image_name, 'signal')
	img = StringIO()
	plt.savefig(img)
	img.seek(0)
	
	return send_file(img, mimetype='image/png')




