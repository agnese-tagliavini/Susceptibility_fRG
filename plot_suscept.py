#!/usr/bin/python
#It performs the plotting for the nonlocal vertex for a fixed FF1 and FF2 (form factors)
#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math
#from agneselibrary.mymath import *

#--------------------------------------SETTINGS ------------------------------------------
#if -1, change overall vertex sign (switch definition)
vert_mul = 1
#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")


fname1 = "dat/carsten_comparison/U2/suscept/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
fname2 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0p3_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_FFT_DIM24_allsymm_withSE.h5"  
#fname2 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta4_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
fname3 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0p3_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname4p5 = "dat/carsten_comparison/U2/suscept/suscept_U2_Beta4p5_PFCB16_HUB_SU2_2D_4PISYMM_maxkpos4_allsymm_withSE.h5"  
#fname5 = "dat/carsten_comparison/U2/suscept/suscept_U2_Beta5_PFCB16_HUB_SU2_2D_4PISYMM_maxkpos4_allsymm_withSE.h5"  
#fname = "dat/" + most_recently_edited
#fname1 = "dat/suscept_U2_Beta1_PFCB128_HUB_SU2_2D_bare.h5"  
#fname2 = "dat/suscept_U2_Beta2_PFCB128_HUB_SU2_2D_bare.h5"  
#fname4 = "dat/suscept_U2_Beta4_PFCB128_HUB_SU2_2D_bare.h5"  
#fname4p5 = "dat/suscept_U2_Beta4p5_PFCB128_HUB_SU2_2D_bare.h5"  
#fname5 = "dat/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_bare.h5"  

if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname3 = fname3.rstrip('\n') # strip newline of fname
f3 = h5py.File(fname3, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals1 = f1["/Params"].attrs.values()
parVals2 = f2["/Params"].attrs.values()
parVals3 = f3["/Params"].attrs.values()

# BETA = 5
UINT1 =  parVals1[0] # follows order in output.cpp
BETA1 =  parVals1[1]
B1 =     parVals1[2]
MU1 =     parVals1[3]

# BETA = 5 with Tprime = 0.3
UINT2 =  parVals2[0] # follows order in output.cpp
BETA2 =  parVals2[1]
B2 =     parVals2[2]
MU2 =     parVals2[3]

# BETA = 5 with Tprime = 0.3 and FFT_DIM = 25
UINT3 =  parVals3[0] # follows order in output.cpp
BETA3 =  parVals3[1]
B3 =     parVals3[2]
MU3 =     parVals3[3]

T_PRIME= 0.3

MAX_KPOS = 4
PATCH_COUNT=4*MAX_KPOS*MAX_KPOS

pi = math.pi

prefact = 1.0
#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
pl.rc('text', usetex=True)
pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"

#-----------------------CONSTRUCT ALL NECESSARY K-PATCHES, ADDING, SUBSTRACTING AND NEG---------------------

#--- functions needed to get from patch to kx/ky/kz-coordinates

def get_xindices( patch ):
    return (patch / (2*MAX_KPOS)) % (2 * MAX_KPOS) - MAX_KPOS +1

def get_yindices( patch ):
    return (patch %(2*MAX_KPOS)) - MAX_KPOS +1

def get_patch(x,y):
    x_patch = (x+MAX_KPOS-1+10000*MAX_KPOS)%(2*MAX_KPOS)
    y_patch = (y+MAX_KPOS-1+10000*MAX_KPOS)%(2*MAX_KPOS)
    return x_patch*2*MAX_KPOS+y_patch

#--- construction of matrices starts
add_k= np.array([[None]*PATCH_COUNT]*PATCH_COUNT)

for i in range(PATCH_COUNT):
    for j in range(PATCH_COUNT):
        xi = get_xindices(i)
        yi = get_yindices(i)
        xj = get_xindices(j)
        yj = get_yindices(j)
        add_k[i,j] = get_patch(xi+xj,yi+yj)


dif_k= np.array([[None]*PATCH_COUNT]*PATCH_COUNT)

for i in range(PATCH_COUNT):
    for j in range(PATCH_COUNT):
        xi = get_xindices(i)
        yi = get_yindices(i)
        xj = get_xindices(j)
        yj = get_yindices(j)
        dif_k[i,j] = get_patch(xi-xj,yi-yj)


inv_k= np.array([None]*PATCH_COUNT)

for i in range(PATCH_COUNT):
        xi = get_xindices(i)
        yi = get_yindices(i)
        inv_k[i] = get_patch(-xi,-yi)
#overwrite the add dif and inv in case of simple patch (only way to get patch_count=16)

def Form_factor(n, k):
    form_factor_arr = np.array([1/math.sqrt(4*pi*pi), 
                                1/math.sqrt(2*pi*pi) * math.cos(pi/MAX_KPOS*get_xindices(k)),
                                1/math.sqrt(2*pi*pi) * math.cos(pi/MAX_KPOS*get_yindices(k)),
                                1/math.sqrt(2*pi*pi) * math.sin(pi/MAX_KPOS*get_xindices(k)),
                                1/math.sqrt(2*pi*pi) * math.sin(pi/MAX_KPOS*get_yindices(k))])
    return form_factor_arr[n] 

Parity = np.array([1,1,1,-1,-1])
#---  Helper functions
def neg( w ):
    return fdim - w - 1

def check_bounds( w1, w2, w1p ):
    if ( w1 < 0 or w1 > fdim - 1 or w2 < 0 or w2 > fdim - 1 or w1p < 0 or w1p > fdim - 1 ):
        return False
    return True

#--------------------------------------Chi PLOTTING ------------------------------------------

print("Plotting suceptibilities ...")


#--- Read

#BETA=5

bosgrid1 = np.array(f1["/suscept_func/bgrid"])
rechi_t1 = vert_mul * np.array(f1["/suscept_func/RE_TRIPLET"])
imchi_t1 = vert_mul * np.array(f1["/suscept_func/IM_TRIPLET"])
rechi_s1 = vert_mul * np.array(f1["/suscept_func/RE_SINGLET"])
imchi_s1 = vert_mul * np.array(f1["/suscept_func/IM_SINGLET"])
rechi_d1 = vert_mul * np.array(f1["/suscept_func/RE_DENSITY"])
imchi_d1 = vert_mul * np.array(f1["/suscept_func/IM_DENSITY"])
rechi_m1 = vert_mul * np.array(f1["/suscept_func/RE_MAGNETIC"])
imchi_m1 = vert_mul * np.array(f1["/suscept_func/IM_MAGNETIC"])

fdim_bos1       = bosgrid1.shape[0]
mom_dim_bos1    = rechi_t1.shape[1]
ffactor_count1  = rechi_t1.shape[2]

#BETA=5 tprime = 0.3

bosgrid2 = np.array(f2["/suscept_func/bgrid"])
rechi_t2 = vert_mul * np.array(f2["/suscept_func/RE_TRIPLET"])
imchi_t2 = vert_mul * np.array(f2["/suscept_func/IM_TRIPLET"])
rechi_s2 = vert_mul * np.array(f2["/suscept_func/RE_SINGLET"])
imchi_s2 = vert_mul * np.array(f2["/suscept_func/IM_SINGLET"])
rechi_d2 = vert_mul * np.array(f2["/suscept_func/RE_DENSITY"])
imchi_d2 = vert_mul * np.array(f2["/suscept_func/IM_DENSITY"])
rechi_m2 = vert_mul * np.array(f2["/suscept_func/RE_MAGNETIC"])
imchi_m2 = vert_mul * np.array(f2["/suscept_func/IM_MAGNETIC"])

fdim_bos2       = bosgrid2.shape[0]
mom_dim_bos2    = rechi_t2.shape[1]
ffactor_count2  = rechi_t2.shape[2]

#BETA=5 tprime = 0.3 FFT_DIM = 24

bosgrid3 = np.array(f3["/suscept_func/bgrid"])
rechi_t3 = vert_mul * np.array(f3["/suscept_func/RE_TRIPLET"]) 
imchi_t3 = vert_mul * np.array(f3["/suscept_func/IM_TRIPLET"]) 
rechi_s3 = vert_mul * np.array(f3["/suscept_func/RE_SINGLET"]) 
imchi_s3 = vert_mul * np.array(f3["/suscept_func/IM_SINGLET"]) 
rechi_d3 = vert_mul * np.array(f3["/suscept_func/RE_DENSITY"]) 
imchi_d3 = vert_mul * np.array(f3["/suscept_func/IM_DENSITY"]) 
rechi_m3 = vert_mul * np.array(f3["/suscept_func/RE_MAGNETIC"])
imchi_m3 = vert_mul * np.array(f3["/suscept_func/IM_MAGNETIC"])

fdim_bos3       = bosgrid3.shape[0]
mom_dim_bos3    = rechi_t3.shape[1]
ffactor_count3  = rechi_t3.shape[2]

#==================================================================================================
#
#       PLOT N1: SUSCEPTIBILITIES IN THE DIFFERENT CHANNELS AS A FUNCTIONS OF TEMPERATURES
#
#==================================================================================================


omega = 0

omega1 = omega + (fdim_bos1 - 1)/2 
omega2 = omega + (fdim_bos2 - 1)/2 
omega3 = omega + (fdim_bos3 - 1)/2 


#q = (0,0)

if (PATCH_COUNT == 16):
    q = 5
elif (PATCH_COUNT == 4):
    q = 0
elif (PATCH_COUNT== 64):
    q = 27

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#SINGLET
rechis_s = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#DENSITY
rechid_s = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#MAGNETIC
rechim_s = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create s*-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_ss = np.array([rechi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_ss = np.array([imchi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#SINGLET              
rechis_ss = np.array([rechi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_ss = np.array([imchi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#DENSITY              
rechid_ss = np.array([rechi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_ss = np.array([imchi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC   

rechim_ss = np.array([rechi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_ss = np.array([imchi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create d-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_d = np.array([rechi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_d = np.array([imchi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#SINGLET 

rechis_d = np.array([rechi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_d = np.array([imchi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#DENSITY  

rechid_d = np.array([rechi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_d = np.array([imchi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC            
rechim_d = np.array([rechi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_d = np.array([imchi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create px-wave array for different temperatures
FF1 = 3
FF2 = 3 

#TRIPLET
rechit_px = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_px = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
#SINGLET              
rechis_px = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_px = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_px = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_px = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_px = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_px = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create py-wave array for different temperatures
FF1 = 4
FF2 = 4 

#TRIPLET
rechit_py = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_py = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_py = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_py = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_py = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_py = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_py = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_py = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


beta_array = np.array([1./BETA3,1./BETA2, 1./BETA1])

#--- Helper functions

def plotchi( use_pl, arrs, arrss, arrd, arrpx, arrpy):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrss,"ro", label = "$s^{*}$")
    pl.plot( beta_array, arrd, "gs", label = "$d_{x^2-y^2}$")
    pl.plot( beta_array, arrpx,"m^", label = "$p_{x}$")
    pl.plot( beta_array, arrpy,"kv", label = "$p_{y}$")
    pl.xlim([0.0, 1.5])
    pl.legend(loc=1, prop={'size':6})
    return



pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_ss, rechit_d, rechit_px, rechit_py) 
pl.ylabel(RE + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_ss, rechis_d, rechis_px, rechis_py) 
pl.ylabel(RE + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_ss, rechid_d, rechid_px, rechid_py) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_ss, rechim_d, rechim_px, rechim_py) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
#pl.tight_layout()

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# RE PLOT FOR POSTER! Density, Magnetic, SC = 1/2(s+t)

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_ss+rechis_ss),0.5*(rechit_d+rechis_d), 0.5*(rechit_px+rechis_px), 0.5*(rechit_py+ rechis_py)) 
pl.ylabel(RE + r"\chi^{SC}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_ss, rechid_d, rechid_px, rechid_py) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_ss, rechim_d, rechim_px, rechim_py) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
#pl.tight_layout()

pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_ss, imchit_d, imchit_px, imchit_py) 
pl.ylabel(IM + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_ss, imchis_d, imchis_px, imchis_py) 
pl.ylabel(IM + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_ss, imchid_d, imchid_px, imchid_py) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_ss, imchim_d, imchim_px, imchim_py) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

def plotchifreq( use_pl, arr1, arr2, arr3):
    #pl.plot( bosgrid1, arr1[:,q, FF1, FF2, 0,0,0,0], "bx", label = "$ beta =5$")
    pl.plot( bosgrid2, arr2[:,q, FF1, FF2, 0,0,0,0], "ro", label = "$FFT_D = 24$")
    pl.plot( bosgrid3, arr3[:,q, FF1, FF2, 0,0,0,0], "gs", label = "$FFT_D = 16$")
    pl.legend(loc=1, prop={'size':6})
    return

FF1 = 0
FF2 = 0

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (0,0)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 1
FF2 = 1

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (0,0)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 3
FF2 = 3

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (\pi,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#====================================================================================
#
#                                   q = (pi,pi)
#
#====================================================================================

if (PATCH_COUNT == 16):
    q = 15
elif (PATCH_COUNT == 4):
    q = 3
elif (PATCH_COUNT== 64):
    q = 63

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#SINGLET             
rechis_s = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#DENSITY             
rechid_s = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#MAGNETIC            
rechim_s = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create p-wave(symmetric) array for different temperatures
FF1 = 3
FF2 = 4 

#TRIPLET
rechit_ps = np.array([rechi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_ps = np.array([imchi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_ps = np.array([rechi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_ps = np.array([imchi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_ps = np.array([rechi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_ps = np.array([imchi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_ps = np.array([rechi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] + rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] + rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] + rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_ps = np.array([imchi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] + imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] + imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] + imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create p-wave(antisymmetric) array for different temperatures
FF1 = 3
FF2 = 4 

#TRIPLET
rechit_pa = np.array([rechi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_pa = np.array([imchi_t3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_pa = np.array([rechi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_pa = np.array([imchi_s3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_pa = np.array([rechi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_pa = np.array([imchi_d3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_pa = np.array([rechi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] - rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] - rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] - rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_pa = np.array([imchi_m3[omega3, q, FF1, FF1, 0, 0, 0, 0] - imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF1, 0, 0, 0, 0] - imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF1, 0, 0, 0, 0] - imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create coskx-wave array for different temperatures
FF1 = 1
FF2 = 1 

#TRIPLET
rechit_cx = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_cx = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                                                              
#SINGLET                                                      
rechis_cx = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_cx = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_cx = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_cx = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                                                              
#MAGNETIC                                                     
rechim_cx = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_cx = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


#Create cosky-wave array for different temperatures
FF1 = 2
FF2 = 2 

#TRIPLET
rechit_cy = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_cy = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_cy = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_cy = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_cy = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_cy = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_cy = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_cy = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])


beta_array = np.array([1./BETA3,1./BETA2, 1./BETA1])

#--- Helper functions

def plotchi( use_pl, arrs, arrcx, arrcy, arrps, arrpa):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrcx,"ro", label = "$cos(kx)$")
    pl.plot( beta_array, arrcy,"gs", label = "$cos(ky)$")
    pl.plot( beta_array, arrps,"m^", label = "$p_{x}+p_y$")
    pl.plot( beta_array, arrpa,"kv", label = "$p_{x}-p_y$")
    pl.xlim([0.0, 1.5])
    pl.legend(loc=1, prop={'size':6})
    return



pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_cx, rechit_cy, rechit_ps, rechit_pa) 
pl.ylabel(RE + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_cx, rechis_cy, rechis_ps, rechis_pa) 
pl.ylabel(RE + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_cx, rechid_cy, rechid_ps, rechid_pa) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_cx, rechim_cy, rechim_ps, rechim_pa) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# RE PLOT POSTER!
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_cx+rechis_cx), 0.5*(rechit_cy+ rechis_cy), 0.5*(rechit_ps+rechis_ps), 0.5*(rechit_pa+rechis_pa)) 
pl.ylabel(RE + r"\chi^{SC}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_cx, rechid_cy, rechid_ps, rechid_pa) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_cx, rechim_cy, rechim_ps, rechim_pa) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#IM

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_cx, imchit_cy, imchit_ps, imchit_pa) 
pl.ylabel(IM + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_cx, imchis_cy, imchis_ps, imchis_pa) 
pl.ylabel(IM + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_cx, imchid_cy, imchid_ps, imchid_pa) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_cx, imchim_cy, imchim_ps, imchim_pa) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


def plotchifreq( use_pl, arr1, arr2, arr3):
    #pl.plot( bosgrid1, arr1[:,q, FF1, FF2, 0,0,0,0], "bx", label = "dat1")
    pl.plot( bosgrid1, arr2[:,q, FF1, FF2, 0,0,0,0], "ro", label = "$FFT_D = 24$")
    pl.plot( bosgrid1, arr3[:,q, FF1, FF2, 0,0,0,0], "gs", label = "$FFT_D = 16$")
    pl.legend(loc=1, prop={'size':6})
    return

FF1 = 0
FF2 = 0

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (\pi,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 1
FF2 = 1

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (\pi,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 3
FF2 = 3

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (\pi,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#===================================================================================================================================
#
#                                       q = (0,pi)
#
#==================================================================================================================================

#q=(0,pi)

if (PATCH_COUNT == 16):
    q = 7
elif (PATCH_COUNT == 4):
    q = 1
elif (PATCH_COUNT== 64):
    q = 31

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#SINGLET             
rechis_s = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#DENSITY             
rechid_s = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                     
#MAGNETIC            
rechim_s = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create cos(kx)-wave array for different temperatures

FF1 = 1
FF2 = 1 

#TRIPLET
rechit_cx = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_cx = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_cx = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_cx = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_cx = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_cx = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_cx = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_cx = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create cos(ky)-wave array for different temperatures

FF1 = 2
FF2 = 2 

#TRIPLET
rechit_cy = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_cy = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLETy             
rechis_cy = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_cy = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_cy = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_cy = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIY             
rechim_cy = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_cy = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create sin(kx)-wave array for different temperatures

FF1 = 3
FF2 = 3 

#TRIPLET
rechit_sx = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_sx = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLET              
rechis_sx = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_sx = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_sx = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_sx = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIC             
rechim_sx = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_sx = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])

#Create sin(ky)-wave array for different temperatures

FF1 = 4
FF2 = 4 

#TRIPLET
rechit_sy = np.array([rechi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchit_sy = np.array([imchi_t3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#SINGLETy             
rechis_sy = np.array([rechi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchis_sy = np.array([imchi_s3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#DENSITY              
rechid_sy = np.array([rechi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchid_sy = np.array([imchi_d3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
                      
#MAGNETIY             
rechim_sy = np.array([rechi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])
imchim_sy = np.array([imchi_m3[omega3, q, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q, FF1, FF2, 0, 0, 0, 0]])



beta_array = np.array([1./BETA3,1./BETA2, 1./BETA1])

#--- Helper functions

def plotchi( use_pl, arrs, arrcx, arrcy, arrsx, arrsy):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrcx,"ro", label = "$\cos(k_x)$")
    pl.plot( beta_array, arrcy,"gs", label = "$\cos(k_y)$")
    pl.plot( beta_array, arrsx,"m^", label = "$\sin(k_x)$")
    pl.plot( beta_array, arrsy,"kv", label = "$\sin(k_y)$")
    pl.xlim([0.0, 1.5])
    pl.legend(loc=1, prop={'size':6})
    return



pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_cx, rechit_cy, rechit_sx, rechit_sy) 
pl.ylabel(RE + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_cx, rechis_cy, rechis_sx, rechis_sy) 
pl.ylabel(RE + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_cx, rechid_cy, rechid_sx, rechid_sy) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_cx, rechim_cy, rechim_sx, rechim_sy) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# PLOT RE POSTER!
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_cx+rechis_cx), 0.5*(rechit_cy+rechis_cy), 0.5*(rechit_sx+rechis_sx), 0.5*(rechit_sy+rechis_sy)) 
pl.ylabel(RE + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_cx, rechid_cy, rechid_sx, rechid_sy) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_cx, rechim_cy, rechim_sx, rechim_sy) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#IM

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_cx, imchit_cy, imchit_sx, imchit_sy) 
pl.ylabel(IM + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_cx, imchis_cy, imchis_sx, imchis_sy) 
pl.ylabel(IM + r"\chi^{s}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_cx, imchid_cy, imchid_sx, imchid_sy) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_cx, imchim_cy, imchim_sx, imchim_sy) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

def plotchifreq( use_pl, arr1, arr2, arr3):
    #pl.plot( bosgrid1, arr1[:,q, FF1, FF2, 0,0,0,0], "bx", label = "$dat1$")
    pl.plot( bosgrid1, arr2[:,q, FF1, FF2, 0,0,0,0], "ro", label = "$FFT_D = 24$")
    pl.plot( bosgrid1, arr3[:,q, FF1, FF2, 0,0,0,0], "gs", label = "$FFT_D = 16$")
    pl.legend(loc=1, prop={'size':6})
    return

FF1 = 0
FF2 = 0

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (0, \pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 1
FF2 = 1

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (0,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

FF1 = 3
FF2 = 3

pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +" $q = (\pi,\pi)$" )

plotchifreq( pl.subplot(2,2,1), rechi_t1, rechi_t2, rechi_t3) 
pl.ylabel(RE + r"\chi^{t}$")
plotchifreq( pl.subplot(2,2,2), rechi_s1, rechi_s2, rechi_s3) 
pl.ylabel(RE + r"\chi^{s}$")
plotchifreq( pl.subplot(2,2,3), rechi_d1, rechi_d2, rechi_d3) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(2,2,4), rechi_m1, rechi_m2, rechi_m3) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_(FF1,FF2)=("+str(FF1)+","+str(FF2)+")_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()



#===============================================================================================================
#
#                                       SELF-ENERGY
#
#===============================================================================================================


fname1 = "dat/carsten_comparison/U2/PREPROC/dat_U2_TP0_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_maxkpos4_allsymm.h5"  
fname2 = "dat/carsten_comparison/U2/PREPROC/dat_U2_TP0p3_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_maxkpos4_allsymm.h5"  
fname3 = "dat/carsten_comparison/U2/PREPROC/dat_U2_TP0p3_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_maxkpos4_allsymm_1FF_FFTDIM=24.h5"  


if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname3 = str(sys.argv[1])

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname1= fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname3 = fname3.rstrip('\n') # strip newline of fname
f3 = h5py.File(fname3, "r")

#--- Read
#BETA=5
siggrid2  = f2["/Sig/fgrid"]
resig2    = prefact * np.array(f2["/Sig/RE"])
imsig2    = prefact * np.array(f2["/Sig/IM"])
mom_dim2  = resig2.shape[1]
fdim2     = resig2.shape[0]
fdimo22   = fdim2/2

#BETA=5
siggrid3  = f3["/Sig/fgrid"]
resig3    = prefact * np.array(f3["/Sig/RE"])
imsig3    = prefact * np.array(f3["/Sig/IM"])
mom_dim3  = resig3.shape[1]
fdim3     = resig3.shape[0]
fdimo23   = fdim3/2

#BETA=1
siggrid1  = f1["/Sig/fgrid"]
resig1    = prefact * np.array(f1["/Sig/RE"])
imsig1    = prefact * np.array(f1["/Sig/IM"])
mom_dim1  = resig1.shape[1]
fdim1     = resig1.shape[0]
fdimo21   = fdim1/2

#==============================================================
#
#       PLOT: SE_REAL AND SE_IMAG FOR DIFFERENT K-POINTS
#
#==============================================================

#--- Helper functions
 #q(0,0)
if PATCH_COUNT==4:
    q = 0
elif PATCH_COUNT==16:
    q = 5
elif PATCH_COUNT==64:
    q = 27

def plotSig( use_pl, arr1, arr2, arr3, string ):
    pl.plot( siggrid1, arr1[:], 'bx',  label="$beta=5$")
    pl.plot( siggrid2, arr2[:], 'go',  label='$beta=5,tp$')
    pl.plot( siggrid3, arr3[:], 'rv',  label='$beta=5,tp$ big')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    pl.legend(loc=2, prop={'size':7})
    return

def plotSigre( use_pl, arr1, arr2, arr3, string ):
    pl.plot( siggrid1, arr1[:], 'bx', label="$beta=5$")
    pl.plot( siggrid2, arr2[:], 'go', label='$beta=5, tp$')
    pl.plot( siggrid3, arr3[:], 'rv', label='$beta=5, tp$ big')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    return

pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (0,0)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig1[:,q,0,0], resig2[:,q,0,0], resig3[:,q,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig1[:,q,0,0], imsig2[:,q,0,0], imsig3[:,q,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--- Helper functions
 #q=(0,pi)

if PATCH_COUNT==4:
    q = 1
elif PATCH_COUNT==16:
    q = 7
elif PATCH_COUNT==64:
    q = 31


pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (0,\pi)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig1[:,q,0,0], resig2[:,q,0,0], resig3[:,q,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig1[:,q,0,0], imsig2[:,q,0,0], imsig3[:,q,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Helper functions
 #q=(pi/2,pi/2)

if PATCH_COUNT==4:
    print "NO (pi/2,pi/2) point for PATCH_COUNT=4"
elif PATCH_COUNT==16:
    q = 10
elif PATCH_COUNT==64:
    q = 45


pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (\pi/2,\pi/2)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig1[:,q,0,0], resig2[:,q,0,0], resig3[:,q,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig1[:,q,0,0], imsig2[:,q,0,0], imsig3[:,q,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(pi2,pi2).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Helper functions
 #q=(pi,pi)

if PATCH_COUNT==4:
    q = 3
elif PATCH_COUNT==16:
    q = 15
elif PATCH_COUNT==64:
    q = 63


pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (\pi,\pi)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig1[:,q,0,0], resig2[:,q,0,0], resig3[:,q,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig1[:,q,0,0], imsig2[:,q,0,0], imsig3[:,q,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()



