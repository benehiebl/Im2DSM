# -*- coding: utf-8 -*-
"""
=================
Im2DSM.py
=================

Version 1.8
Author: Bene Hiebl
Date: 11.10.2022

Monoplotting tool for projection of single channel raster images to DSM.
Using only first 3 GCPs in list!

Input requirements:
•	DSM
•	Image (.asc, .tif)
•	Camera Location 
•	Camera Parameters (focal length, sensor size or FOV)
•	Ground Control Points or precalculated pitch, yaw, roll angles

Define Input in Im2DSM_Parameters.txt


"""
import os

import matplotlib.pyplot as plt

import numpy as np
import pandas
import math as m

from osgeo import gdal
import rasterio

import warnings
warnings.filterwarnings('ignore')

#############################################################################################
### FUNCTIONS ###
##############################################################################################

def file_list(dir,file_format):
    "creates string list of files in specific folder (format: foldername, no /)"
    
    image_dir = os.scandir(dir+"/")
    imagenames =[]
    for file in image_dir:
        name = file.name
        if name.endswith(file_format):
            imagenames.append(name)
    return imagenames

def read_parameterfile(dir_path, file_name):
    "Reads Parameterfile from txt and give information to dictionary, returns dictionary"
    
    pfile = open(dir_path+"/"+file_name,"r")

    # create dictionary from parameters
    try:
        pdict = {}
        for line in pfile:
    
            # delete # comments 
            if line.startswith("#") == False:
                # delete certain spaces that are not part of the file
                line = line.strip(" ")
                # delete \n
                line = line.strip("\n")
                # tabulator as seperator 
                sep_line = line.split("\t")
                sep_line = list(sep_line)
                # values to dictionary
                pdict[sep_line[0]] = sep_line[2]
        
        pfile.close()
    except:
        print("Something's wrong with the Parameter file! Probably you missed a /t (Tab)...")
        
    # convert numbers to float    
    for i in pdict:
        try:
            pdict[i] = float(pdict[i])
        except:
            pdict[i] = str(pdict[i])
    
        
    return pdict

def out_directory (directory,new_folder):
    "creates output directory in same folder as .py file"
    
    if new_folder.startswith("/"):
        out_dir = directory + new_folder
    else:
        out_dir = directory + "/" + new_folder
        
    try:
        os.mkdir(out_dir)
    except OSError:
        print ("Creation of the directory /"+new_folder+" failed")
    else:
        print ("Successfully created the directory /"+new_folder)
    return out_dir

def pyplot_arrays (data,v_min,v_max,title):
    "creates Simple mesh plots"
    
    cf = plt.imshow(data,vmin=v_min,vmax=v_max,extent = [0,np.shape(data)[1],0,np.shape(data)[0]])
    plt.colorbar(cf)
    plt.title(title)
    plt.show()
    
def hill_slope_aspect (out_dir, dsm_name):
    "creates slope, aspect, hillshade as Gtiff"
    
    gdal.DEMProcessing(out_dir+"/slope.tif", dsm_name, "slope")
    gdal.DEMProcessing(out_dir+"/aspect.tif", dsm_name, "aspect")
    gdal.DEMProcessing(out_dir+"/hillshade.tif", dsm_name, "hillshade")

### convert array to Gtiff
def array_to_Gtiff(out_name,data_array,out_dir,geotransform,projection):
    "Saves np.array to Geotiff, note: out_dir has to be defined previously"
    
    drv = gdal.GetDriverByName("GTiff")
    width = np.shape(data_array)[1]
    height = np.shape(data_array)[0]
    ds = drv.Create(out_dir+"/"+out_name+".tif", width, height, 1, gdal.GDT_Float32)
    ds.SetGeoTransform(geotransform)
    ds.SetProjection(projection)
    ds.GetRasterBand(1).WriteArray(data_array)

### get raster coordinates from pixel in array of certain geotransformation
def rast_coordinates(i,j,dsm,dsm_array):
    "computes raster coordinates from array coordinates, returns list=(x,y,z)"
    
    x = i*dsm.GetGeoTransform()[1] + dsm.GetGeoTransform()[0]
    y = j*dsm.GetGeoTransform()[5] + dsm.GetGeoTransform()[3]
    z = dsm_array[j,i]
    loc = (x, y, z)
    return loc

### 
def pix_coordinates(x,y,dsm,dsm_array):
    "computes array coordinates from raster coordinates, returns list=(i,j,z)"
    
    i = (x - dsm.GetGeoTransform()[0]) / dsm.GetGeoTransform()[1]
    j = (y - dsm.GetGeoTransform()[3]) / dsm.GetGeoTransform()[5]
    z = dsm_array[round(j),round(i)]
    loc = (i,j,z)
    return loc

###
def coordinate_angles(cam_location, point_location):
    "creates horizontal and vertical angle between points, returns list of angles"
    
    # vertical angle
    ground_dist = m.sqrt((cam_location[0]-point_location[0])**2 + (cam_location[1]-point_location[1])**2)
    height_diff = point_location[2] - cam_location[2]
    h_degree = m.atan(height_diff/ground_dist) #* (180/m.pi)
    h_degree = m.degrees(h_degree)
    # horizontal angle (in coordinate system)
    x_degree = m.atan((point_location[1]-cam_location[1])/(point_location[0]-cam_location[0])) #* (180/m.pi)
    x_degree = m.degrees(x_degree)
    angles = (x_degree,h_degree)
    return angles


def viewshed_array (out_name,cam_location, dsm_name):
    "creates viewshed mask as array from dsm and cam location and saves to OUTPUT location, returns viewshed array"
    
    os.system("C:/Users/beneh/Anaconda3/envs/rs/Library/bin/gdal_viewshed.exe -ox " + str(cam_location[0]) + " -oy " + str(cam_location[1]) + " -oz " + "1" + " -vv 1 " + dsm_name + " " + out_name)      
    #vis_mask = gdal.Open(out_name+".tif")
    vis_mask = gdal.Open(out_name)
    viewshed = np.array(vis_mask.GetRasterBand(1).ReadAsArray())   
    return viewshed


##############################################################################################
### Input ###
#############################################################################################

### Read Parameter file and give to dictionary
dir_path = os.path.dirname(os.path.realpath("Im2DSM.py"))
pdict = read_parameterfile(dir_path, "Im2DSM_Parameters.txt")

# working directory
in_dir = pdict["IN_DIRECTORY"]

# output directory
out_dir = pdict["OUT_DIRECTORY"]

# filenames
dsm_name = pdict["DSM_INPUT"]
viewshed_name = pdict["VIEWSHED"]
asc_skip = pdict["ASC_SKIP"]
file_format = pdict["FILE_FORMAT"]

# Camera location coordinates and specific camera parameters (fov = (angle horizontal, angle vertical))
cam = (pdict["CAMLOCATION_X"],pdict["CAMLOCATION_Y"],pdict["CAMLOCATION_Z"])
focal = pdict["FOCAL_LENGTH"]
sensor_x = pdict["SENSOR_X"]
sensor_y = pdict["SENSOR_Y"]
fov_x = pdict["FOV_SENSOR_X"]
fov_y = pdict["FOV_SENSOR_Y"]

#fov = (pdict["FOV_CAM_X"],pdict["FOV_CAM_Y"])
distortion_a = pdict["DISTORTION"]
dist_zero = pdict["DIST_ZERO"]
dist_mode = pdict["DIST_MODE"]
adjust_left = pdict["ADJUST_LEFT"]
adjust_up = pdict["ADJUST_UP"]
adjust_roll = pdict["ADJUST_ROLL"]
colin = pdict["CO_LINEAR"]

# INS
INS_angle_up = pdict["ANGLE_UP"]
INS_angle_left = pdict["ANGLE_LEFT"]
INS_angle_rot = pdict["ANGLE_ROT"]

# Reference points for camera orientation (if no INS data is available)
gcp_file = pdict["GCP_FILE"]
accuracy_mode = pdict["ACCURACY"]


print("reading parameter file.")
print(file_format," ",focal," ",sensor_x," ",sensor_y," ",distortion_a," ",adjust_left," ",adjust_up," ",adjust_roll)


###########################################################################################
### Start Script #######
###########################################################################################
print ("Reading data.")

# create output folder
foldername = in_dir.split("/")
out_dir = out_directory(out_dir,foldername[-1])

### Get names of image input files
imagenames = file_list(in_dir, file_format)

### read data DSM and Pic ###
if imagenames[0].endswith(".asc"):
    file = pandas.read_table(in_dir+"/"+imagenames[0],delimiter="\t", decimal = ".", skiprows=int(asc_skip), header=None)    
    th_ar = np.array(file)
elif imagenames[0].endswith(".tif"):
    with rasterio.open(in_dir+"/"+imagenames[0],"r") as img:
        th_ar = img.read(1)
    # with gdal.Open(in_dir+"/"+imagenames[0]) as img:
    #     th_ar = np.array(img.GetRasterBand(1).ReadAsArray())

# read dsm
dsm = gdal.Open(dsm_name)

### create arrays ###
dsm_ar = np.array(dsm.GetRasterBand(1).ReadAsArray())

# read GCP file
gcp = pandas.read_csv(gcp_file)
# add dsm heihgt col to gcp s
dsm_height = []
for i in range(0,len(gcp["Name"])):
    coord = pix_coordinates(gcp["X"][i], gcp["Y"][i], dsm, dsm_ar)
    dsm_height.append(coord[2])                              
gcp["dsm_height"] = dsm_height   

# picture size
pic_y = th_ar.shape[0]
pic_x = th_ar.shape[1]

#####################################################################################################
### Getting orientation of camera by getting angles of upper left pixel and function for pix deviation ###
#####################################################################################################
### camera orientation only for straight angles
print ("set camera orientation parameters.")

# calculate fov of camera
if fov_x != "NA":
    fov = fov_x, fov_y
elif fov_x == "NA":
    if focal != "NA":
        fov = np.degrees(2 * np.arctan((sensor_x/2) / focal)), np.degrees(2 * np.arctan((sensor_y/2) / focal))
        #fov = fov1[0] * (1 + distortion_a*(sensor_x/sensor_y)**3 + distortion_a*1), fov1[1] * (1 + distortion_a*1**3 + distortion_a*1)
    
    elif focal == "NA":
        ### FOV from GCP
        base_vector = [gcp["X"][0]-cam[0], gcp["Y"][0]-cam[1], gcp["dsm_height"][0]-cam[2]]
        x_fov = []
        y_fov = []
        for nr in range(1,len(gcp)):
            sec_vector = [gcp["X"][nr]-cam[0], gcp["Y"][nr]-cam[1], gcp["dsm_height"][nr]-cam[2]]
            dot = base_vector[0]*sec_vector[0] + base_vector[1]*sec_vector[1] + base_vector[2]*sec_vector[2]
            lbase = np.sqrt(base_vector[0]**2 + base_vector[1]**2 + base_vector[2]**2)
            lsec = np.sqrt(sec_vector[0]**2 + sec_vector[1]**2 + sec_vector[2]**2)
            alpha = np.arccos(dot / (lbase*lsec))
            
            x = gcp["pic_X"][nr]-gcp["pic_X"][0]
            y = gcp["pic_Y"][nr]-gcp["pic_Y"][0]
            
            # fov in x direction
            dx = np.sqrt(x**2 + y**2)
            dxf = np.sqrt(((y / x) * pic_x)**2 + pic_x**2) 
            
            alpha_f = (dxf / dx) * alpha
            
            fov_x = np.sqrt(((x / dx) * alpha_f)**2)
            x_fov.append(fov_x)
            
            #fov in y direction
            dyf = np.sqrt(((x / y) * pic_y)**2 + pic_y**2)
            alphay_f = (dyf / dx) * alpha
            
            fov_y = np.sqrt(((y / dx) * alphay_f)**2)
            y_fov.append(fov_y)
            
        fov = np.degrees(np.nanmean(x_fov)), np.degrees(np.nanmean(y_fov))
        #calculate new fov
        #fov = fov1[0] * (1 + distortion_a*(sensor_x/sensor_y)**3 + distortion_a*1), fov1[1] * (1 + distortion_a*1**3 + distortion_a*1)


print("Used FOV: " + str(fov))
    

# Calculate camera height from DOM
if viewshed_name == "None":
    cam_array_loc = pix_coordinates(cam[0], cam[1], dsm, dsm_ar)
    cam = list(cam)
    cam[2] = cam_array_loc[2]

# coordinates of image center
center_coordinates = [th_ar.shape[1]/2,th_ar.shape[0]/2]

# camera orientation

if INS_angle_up == "NA":
    
    # calculate distance camera to image plane in px
    d = (th_ar.shape[1]/2) / np.tan(np.radians(fov[0]) / 2)
    
    angles_left = []
    angles_up = []
    angles_roll = []
    
    ###  calculate viewing direction, pitch and yaw
    
    # distance camera to GCP on image plane in px for 3 GCPs
    img1 = np.sqrt((np.sqrt((list(gcp["pic_X"])[0]-th_ar.shape[1]/2)**2))**2 + d**2 + (np.sqrt((list(gcp["pic_Y"])[0]-th_ar.shape[0]/2)**2))**2)
    img2 = np.sqrt((np.sqrt((list(gcp["pic_X"])[1]-th_ar.shape[1]/2)**2))**2 + d**2 + (np.sqrt((list(gcp["pic_Y"])[1]-th_ar.shape[0]/2)**2))**2)
    img3 = np.sqrt((np.sqrt((list(gcp["pic_X"])[2]-th_ar.shape[1]/2)**2))**2 + d**2 + (np.sqrt((list(gcp["pic_Y"])[2]-th_ar.shape[0]/2)**2))**2)
    
    # calculate real distances in WGS84 to image plane
    gcp1 = [list(gcp["X"])[0]-cam[0],list(gcp["Y"])[0]-cam[1],list(gcp["dsm_height"])[0]-cam[2]] 
    gcp1_new = (gcp1 / np.linalg.norm(gcp1)) * img1
    gcp2 = [list(gcp["X"])[1]-cam[0],list(gcp["Y"])[1]-cam[1],list(gcp["dsm_height"])[1]-cam[2]] 
    gcp2_new = (gcp2 / np.linalg.norm(gcp2)) * img2
    gcp3 = [list(gcp["X"])[2]-cam[0],list(gcp["Y"])[2]-cam[1],list(gcp["dsm_height"])[2]-cam[2]] 
    gcp3_new = (gcp3 / np.linalg.norm(gcp3)) * img3
    
    imggcp1 = [list(gcp["pic_X"])[0]-th_ar.shape[1]/2, list(gcp["pic_Y"])[0]-th_ar.shape[0]/2, d]
    imggcp3 = [list(gcp["pic_X"])[2]-th_ar.shape[1]/2, list(gcp["pic_Y"])[2]-th_ar.shape[0]/2, d]
    
    img_angle = np.arccos(np.dot(imggcp1,imggcp3) / (np.sqrt(np.dot(imggcp1,imggcp1)) * np.sqrt(np.dot(imggcp3,imggcp3))))
    
    
    real_angle = np.arccos(np.dot(gcp1,gcp3) / (np.sqrt(np.dot(gcp1,gcp1)) * np.sqrt(np.dot(gcp3,gcp3))))
    
    # calculate vectors between 3 GCPs to get location of image plane in WGS84
    a = np.array(gcp2_new) - np.array(gcp1_new)
    b = np.array(gcp3_new) - np.array(gcp1_new)
    
    # calculate viewing direction as normal vector on image plane
    z = np.cross(a,b)
    z_angle = np.arccos(np.dot(z,gcp1_new) / (np.sqrt(np.dot(z,z)) * np.sqrt(np.dot(gcp1_new,gcp1_new))))
    if z_angle > np.radians(90):
        z = np.cross(b,a)
    
    #length adapted z vector (image axes, viewing direction)
    z_norm = (z / np.linalg.norm(z)) * d
    
    #calculate pitch and yaw from z            
    angle_left = np.degrees(np.arctan(z_norm[1] / z_norm[0]))
    angles_left.append(angle_left)
    angle_up = np.degrees(np.arccos(np.sqrt(z_norm[0]**2 + z_norm[1]**2) / np.sqrt(z_norm[0]**2 + z_norm[1]**2 + z_norm[2]**2)))
    if z_norm[2] < 0:
        angles_up.append(-angle_up)
    else:
        angles_up.append(angle_up)
    
    
    #### implement roll angle
    
    # image center to GCP vector
    img_a = [list(gcp["pic_X"])[0] - th_ar.shape[1]/2, list(gcp["pic_Y"])[0] - th_ar.shape[0]/2]
    img_b = [list(gcp["pic_X"])[1] - th_ar.shape[1]/2, list(gcp["pic_Y"])[1] - th_ar.shape[0]/2]
    img_c = [list(gcp["pic_X"])[2] - th_ar.shape[1]/2, list(gcp["pic_Y"])[2] - th_ar.shape[0]/2]
    
    # image x  axis vector (pointing towards increasing x (right))
    img_x = [1,0]
    # angle between x axis vector and center to GCP vector
    img_ax_angle = np.arccos(np.dot(img_a,img_x) / (np.sqrt(np.dot(img_a,img_a)) * np.sqrt(np.dot(img_x,img_x))))
    img_bx_angle = np.arccos(np.dot(img_b,img_x) / (np.sqrt(np.dot(img_b,img_b)) * np.sqrt(np.dot(img_x,img_x))))
    img_cx_angle = np.arccos(np.dot(img_c,img_x) / (np.sqrt(np.dot(img_c,img_c)) * np.sqrt(np.dot(img_x,img_x))))
    
    # image center to GCP vector in WGS84 on projected image plane
    wgs_a = np.array(gcp1_new) - z_norm
    wgs_b = np.array(gcp2_new) - z_norm
    wgs_c = np.array(gcp3_new) - z_norm
    
    
    # image plane x axis direction in WGS84 (pointing towards increasing x in image(right))
    if z[0] < 0:
        x_new = [z_norm[1] / -z_norm[0], 1, 0]
    else:
        x_new = [z_norm[1] / z_norm[0], -1, 0]

    # angle between x_new (image x axis in WGS84) and center to GCP vector
    az_angle = np.arccos(np.dot(wgs_a,x_new) / (np.sqrt(np.dot(wgs_a,wgs_a)) * np.sqrt(np.dot(x_new,x_new))))
    bz_angle = np.arccos(np.dot(wgs_b,x_new) / (np.sqrt(np.dot(wgs_b,wgs_b)) * np.sqrt(np.dot(x_new,x_new)))) 
    cz_angle = np.arccos(np.dot(wgs_c,x_new) / (np.sqrt(np.dot(wgs_c,wgs_c)) * np.sqrt(np.dot(x_new,x_new)))) 
    
    # angle roll as difference betweenangle in image and angle in reality
    if img_a[1] > 0:            
        angle_roll_a = np.degrees(az_angle - img_ax_angle)
    else:
        angle_roll_a = np.degrees(img_ax_angle - az_angle)
    
    if img_b[1] > 0:
        angle_roll_b = np.degrees(bz_angle - img_bx_angle)
    else:
        angle_roll_b = np.degrees(img_bx_angle - bz_angle)
    
    if img_c[1] > 0:
        angle_roll_c = np.degrees(cz_angle - img_cx_angle)
    else:
        angle_roll_c = np.degrees(img_cx_angle - cz_angle)
    
    angle_roll = (angle_roll_a + angle_roll_b + angle_roll_c) / 3
    
    
    
    ### calculate final angles as mean of all GCPs and add manual adjustments
    angle_left = np.nanmean(angles_left) + fov[0]/2 + adjust_left
    angle_up = np.nanmean(angles_up) + fov[1]/2 + adjust_up   
    angle_roll = angle_roll + adjust_roll
    
    
else:
    angle_up = INS_angle_up
    angle_left = INS_angle_left
    angle_rot = INS_angle_rot
    
    print("   used INS to calculate camera orientation")

print ("   yaw =   " + str(angle_left - fov[0]/2) + "\n   pitch =   " + str(angle_up - fov[1]/2) + "\n   roll =   " + str(angle_roll))




###################################################################################################
### filling new raster with thermo values
###################################################################################################
print ("Processing image parameters to raster.")

# open DSM raster and convert to array
dsm_th = dsm 
dsm_th_ar = np.array(dsm_th.GetRasterBand(1).ReadAsArray())


## Calculate XYZ arrays
transform = dsm.GetGeoTransform()
y_ar = np.zeros(np.shape(dsm_th_ar))          
for i in range(np.shape(dsm_th_ar)[0]):
    y_ar[i,:] = i*transform[5] + transform[3]

x_ar = np.zeros(np.shape(dsm_th_ar))
for j in range(np.shape(dsm_th_ar)[1]):
    x_ar[:,j] = j*transform[1] + transform[0]

## Pixelkoordinate from pitch and yaw angle
p_y = np.round((angle_up - (np.degrees(np.arctan((dsm_th_ar - cam[2]) / (np.sqrt((cam[0]-x_ar)**2 + (cam[1]-y_ar)**2)))))) * (pic_y/fov[1]))
p_x = np.round((angle_left - (np.degrees(np.arctan((y_ar-cam[1])/(x_ar-cam[0]))))) * (pic_x/fov[0]))

### correct with roll angle
r = np.sqrt((center_coordinates[0]-p_x)**2 + (center_coordinates[1]-p_y)**2)
# calculate new angle
alpha = np.where(p_x == center_coordinates[0], 0, (np.arctan((p_y-center_coordinates[1])/(p_x-center_coordinates[0]))) + (angle_roll/(180/np.pi))) 

# calculate image arrays withoout distortion
th_ax_new = np.where(p_x > center_coordinates[0], np.round((r*np.cos(alpha)) + center_coordinates[0]), center_coordinates[0] - np.round(r*np.cos(alpha)))
th_ay_new = np.where(p_x > center_coordinates[0], np.round((r*np.sin(alpha)) + center_coordinates[1]), center_coordinates[1] - np.round(r*np.sin(alpha)) )

# th_ax_new = np.where(alpha == 0, 0 , th_ax_new)
# test = np.where(alpha == 0, th_ay_new, np.nan)
# plt.imshow(test)
# plt.show()

### calculate distortion
if distortion_a != 0:
    # calculate radius from image center
    r_ar = np.where(th_ax_new != -99999,np.sqrt((th_ax_new-center_coordinates[0])**2 + (th_ay_new-center_coordinates[1])**2), th_ax_new)
    norm_r = np.where(th_ax_new != -99999, np.sqrt(r_ar**2), th_ax_new)
    
    # calculate new image coordinates in real world
    dist_zero = (pic_y/2) * dist_zero
    
    if dist_mode == "SIN":
        dist_x = distortion_a*np.sin(norm_r*(np.pi**2/(dist_zero)))
    elif dist_mode == "SIN_SQRT":
        dist_x = distortion_a*np.sin(np.sqrt(norm_r*(np.pi**2/(dist_zero))))
    else:
        print("distortion correction was set, but wrong mode was chosen! (Use SIN or SIN_SQRT)")
        exit()
    
    plt.imshow(dist_x)
    plt.show()
    
    th_ax = np.where(th_ax_new <= center_coordinates[0], th_ax_new - (np.sqrt((th_ax_new - center_coordinates[0])**2)*dist_x), th_ax_new + (np.sqrt((th_ax_new - center_coordinates[0])**2)*dist_x))
    th_ay = np.where(th_ay_new <= center_coordinates[1], th_ay_new - (np.sqrt((th_ay_new - center_coordinates[1])**2)*dist_x), th_ay_new + (np.sqrt((th_ay_new - center_coordinates[1])**2)*dist_x))
    
    th_ax = np.where(alpha == 0, -99999, th_ax)
    
elif distortion_a == 0:
    th_ax = th_ax_new
    th_ay = th_ay_new

th_ax = np.nan_to_num(th_ax,nan=-99999)                     # delete np.nan
th_ay = np.nan_to_num(th_ay,nan=-99999)





print ("   done.")
### apply viewshed mask to pixelvalues ######################################################################################
print ("mask with viewshed.")

if viewshed_name == "None":
    viewshed_mask = viewshed_array(out_dir+"/viewshed.tif",cam, dsm_name)
else:
    viewshed_tif = gdal.Open(viewshed_name)
    viewshed_mask = np.array(viewshed_tif.GetRasterBand(1).ReadAsArray())

# apply mask to x and y arrays
th_ax = np.where(viewshed_mask == 1, th_ax, -99999).astype(int)
th_ay = np.where(viewshed_mask == 1, th_ay, -99999).astype(int)

plt.imshow(th_ax)
plt.show()

#####################################################################################################################
#### ACCURACY  ########################No documentation availabl##########################
if accuracy_mode != "NA":
    
    print("Calculate accuracy for GCPs.")
    ## check accuracy of GCPs
    check_ar = [th_ax,th_ay]
    check_ar = np.array(check_ar)
    
    check_gcp = []
    for nr in range(0,len(gcp["Name"])):
        val_pic = [pic_y-int(gcp["pic_Y"][nr]),int(gcp["pic_X"][nr])]
        check_gcp.append(val_pic)
        
    
    
    accuracy = []
    
    for i in range(0,check_ar.shape[1]):
        for j in range(0,check_ar.shape[2]):
           
            values = [int(check_ar[1,i,j]),int(check_ar[0,i,j])] 
                  
            if values in check_gcp:                                                           # check if in GCP list
                             
                loc = rast_coordinates(j,i, dsm_th, th_ax)           
                                                      # get position of projected GCP               
                count = np.nan
                for nr in range(0, len(check_gcp)):                                               # get position of value
                    if values == check_gcp[nr]:
                        count = nr
                                         
               
                          
                dist = np.sqrt((loc[0] - gcp["X"][count])**2 + (loc[1] - gcp["Y"][count])**2)
                vector_x, vector_y = loc[0]-gcp["X"][count], loc[1]-gcp["Y"][count]
                
                vals = [gcp["Name"][count],dist,loc[0],loc[1],values[0],values[1],i,j,vector_x,vector_y]
                
                accuracy.append(vals)
                
    accuracy = pandas.DataFrame(accuracy,columns=["Name","dist","N","E","x","y","i","j","Vector_x","Vector_y"]) 
    accuracy["dist"] = accuracy["dist"].astype("float64")         
    
    # save to csv
    accuracy.to_csv(out_dir+"/accuracy.csv", sep=",")
    
    # plot accuracy of GCPs
    plt.plot(accuracy["Name"],accuracy["dist"],"o") 
    plt.show()  

#######################################################################################################################
###### PROJECTING  ########################################################
             
###create raster for every image
print("Start projecting images.")

for name in imagenames:
    
    if name.endswith(".asc"):
        file = pandas.read_table(in_dir+"/"+name,delimiter="\t", decimal = "," , skiprows=int(asc_skip), header=None)    
        im_ar = np.array(file)
    elif name.endswith(".tif"):
        with rasterio.open(in_dir+"/"+name,"r") as img:
            im_ar = img.read(1)
            
    dsm_im = np.zeros(dsm_th_ar.shape)
    for j in range(dsm_th.RasterYSize):
        for i in range(dsm_th.RasterXSize):
            # get pixelcoordinates of image
            # x = int(th_ax[j,i])
            # y = int(th_ay[j,i])
            # give pixelvalue to array
            if 0 < th_ay[j,i] < pic_y and 0 < th_ax[j,i] < pic_x:
                dsm_im[j,i] = im_ar[th_ay[j,i],th_ax[j,i]]
            else:
                dsm_im[j,i] = -99999
    
    #dsm_im = np.where(0 < th_ay < pic_y and 0 < th_ax < pic_x, im_ar[th_ay,th_ax], -99999)
    
    ## Save image array as gtiff    
    array_to_Gtiff(name[:-4]+"_raster", dsm_im, out_dir, dsm.GetGeoTransform(), dsm.GetProjection())
    
    print("   "+name+" done.")
    



print ("projection done.")


