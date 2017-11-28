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
#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

most_recently_edited = run("ls -Art dat/ | tail -n 1")


fname1 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta1_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
fname2 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta2_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
fname4 = "dat/carsten_comparison/U2/suscept/OMFL/suscept_U2_TP0_Beta4_PFCB128_HUB_SU2_2D_4PISYMM_mkp2_ALLSYMM_OMFL.h5"  
fname4p5 = "dat/carsten_comparison/U2/suscept/OMFL/suscept_U2_TP0_Beta4_PFCB128_HUB_SU2_2D_4PISYMM_mkp2_ALLSYMM_OMFL.h5"  
fname5 = "dat/carsten_comparison/U2/suscept/OMFL/suscept_U2_TP0_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp2_ALLSYMM_OMFL.h5"  

#fname1 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta1_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname2 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta2_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname4 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta4_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname4p5 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta4p5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname5 = "dat/carsten_comparison/U2/suscept/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname = "dat/" + most_recently_edited
#fname1 = "dat/suscept_U2_Beta1_PFCB128_HUB_SU2_2D_bare.h5"  
#fname2 = "dat/suscept_U2_Beta2_PFCB128_HUB_SU2_2D_bare.h5"  
#fname4 = "dat/suscept_U2_Beta4_PFCB128_HUB_SU2_2D_bare.h5"  
#fname4p5 = "dat/suscept_U2_Beta4p5_PFCB128_HUB_SU2_2D_bare.h5"  
#fname5 = "dat/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_bare.h5"  
#suscept_U2_Beta1_PFCB128_HUB_SU2_2D_4PISYMM_maxkpos4_allsymm_withoutSE.h5
if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname4 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname4p5 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname5 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname4 = fname4.rstrip('\n') # strip newline of fname
f4 = h5py.File(fname4, "r")

fname4p5 = fname4p5.rstrip('\n') # strip newline of fname
f4p5 = h5py.File(fname4p5, "r")

fname5 = fname5.rstrip('\n') # strip newline of fname
f5 = h5py.File(fname5, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals1 = f1["/Params"].attrs.values()
parVals2 = f2["/Params"].attrs.values()
parVals4 = f4["/Params"].attrs.values()
parVals4p5 = f4p5["/Params"].attrs.values()
parVals5 = f5["/Params"].attrs.values()

# BETA = 1
UINT1 =  parVals1[0] # follows order in output.cpp
BETA1 =  parVals1[1]
B1 =     parVals1[2]
MU1 =     parVals1[3]

# BETA = 2
UINT2 =  parVals2[0] # follows order in output.cpp
BETA2 =  parVals2[1]
B2 =     parVals2[2]
MU2 =     parVals2[3]

# BETA = 4
UINT4 =  parVals4[0] # follows order in output.cpp
BETA4 =  parVals4[1]
B4 =     parVals4[2]
MU4 =     parVals4[3]

# BETA = 4.5
UINT4p5 =  parVals4p5[0] # follows order in output.cpp
BETA4p5 =  parVals4p5[1]
B4p5 =     parVals4p5[2]
MU4p5 =     parVals4p5[3]

# BETA = 5
UINT5 =  parVals5[0] # follows order in output.cpp
BETA5 =  parVals5[1]
B5 =     parVals5[2]
MU5 =    parVals5[3]

T_PRIME=  0 #parVals[4]

MAX_KPOS1 = 2
PATCH_COUNT1=4*MAX_KPOS1*MAX_KPOS1

MAX_KPOS2 = 4
PATCH_COUNT2=4*MAX_KPOS2*MAX_KPOS2

pi = math.pi

#vert_mul = 1.0/4.0/pi/pi
vert_mul = 1.0
prefact = 1.0
#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=10) 
pl.rc('ytick', labelsize=10) 
pl.rc('text', usetex=True)
pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

#-----------------------CONSTRUCT ALL NECESSARY K-PATCHES, ADDING, SUBSTRACTING AND NEG---------------------

#--- functions needed to get from patch to kx/ky/kz-coordinates

def get_xindices( patch ):
    return (patch / (2*MAX_KPOS1)) % (2 * MAX_KPOS1) - MAX_KPOS1 +1

def get_yindices( patch ):
    return (patch %(2*MAX_KPOS1)) - MAX_KPOS1 +1

def get_patch(x,y):
    x_patch = (x+MAX_KPOS1-1+10000*MAX_KPOS1)%(2*MAX_KPOS1)
    y_patch = (y+MAX_KPOS1-1+10000*MAX_KPOS1)%(2*MAX_KPOS1)
    return x_patch*2*MAX_KPOS1+y_patch

#--- construction of matrices starts
add_k= np.array([[None]*PATCH_COUNT1]*PATCH_COUNT1)

for i in range(PATCH_COUNT1):
    for j in range(PATCH_COUNT1):
        xi = get_xindices(i)
        yi = get_yindices(i)
        xj = get_xindices(j)
        yj = get_yindices(j)
        add_k[i,j] = get_patch(xi+xj,yi+yj)


dif_k= np.array([[None]*PATCH_COUNT1]*PATCH_COUNT1)

for i in range(PATCH_COUNT1):
    for j in range(PATCH_COUNT1):
        xi = get_xindices(i)
        yi = get_yindices(i)
        xj = get_xindices(j)
        yj = get_yindices(j)
        dif_k[i,j] = get_patch(xi-xj,yi-yj)


inv_k= np.array([None]*PATCH_COUNT1)

for i in range(PATCH_COUNT1):
        xi = get_xindices(i)
        yi = get_yindices(i)
        inv_k[i] = get_patch(-xi,-yi)
#overwrite the add dif and inv in case of simple patch (only way to get patch_count=16)

def Form_factor(n, k):
    form_factor_arr = np.array([1/math.sqrt(4*pi*pi), 
                                1/math.sqrt(2*pi*pi) * math.cos(pi/MAX_KPOS1*get_xindices(k)),
                                1/math.sqrt(2*pi*pi) * math.cos(pi/MAX_KPOS1*get_yindices(k)),
                                1/math.sqrt(2*pi*pi) * math.sin(pi/MAX_KPOS1*get_xindices(k)),
                                1/math.sqrt(2*pi*pi) * math.sin(pi/MAX_KPOS1*get_yindices(k))])
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

#BETA=1

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

#BETA=2

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

#BETA=4

bosgrid4 = np.array(f4["/suscept_func/bgrid"])
rechi_t4 = vert_mul * np.array(f4["/suscept_func/RE_TRIPLET"])
imchi_t4 = vert_mul * np.array(f4["/suscept_func/IM_TRIPLET"])
rechi_s4 = vert_mul * np.array(f4["/suscept_func/RE_SINGLET"])
imchi_s4 = vert_mul * np.array(f4["/suscept_func/IM_SINGLET"])
rechi_d4 = vert_mul * np.array(f4["/suscept_func/RE_DENSITY"])
imchi_d4 = vert_mul * np.array(f4["/suscept_func/IM_DENSITY"])
rechi_m4 = vert_mul * np.array(f4["/suscept_func/RE_MAGNETIC"])
imchi_m4 = vert_mul * np.array(f4["/suscept_func/IM_MAGNETIC"])

fdim_bos4       = bosgrid4.shape[0]
mom_dim_bos4    = rechi_t4.shape[1]
ffactor_count4  = rechi_t4.shape[2]

#BETA=4p5

bosgrid4p5 = np.array(f4p5["/suscept_func/bgrid"])
rechi_t4p5 = vert_mul * np.array(f4p5["/suscept_func/RE_TRIPLET"])
imchi_t4p5 = vert_mul * np.array(f4p5["/suscept_func/IM_TRIPLET"])
rechi_s4p5 = vert_mul * np.array(f4p5["/suscept_func/RE_SINGLET"])
imchi_s4p5 = vert_mul * np.array(f4p5["/suscept_func/IM_SINGLET"])
rechi_d4p5 = vert_mul * np.array(f4p5["/suscept_func/RE_DENSITY"])
imchi_d4p5 = vert_mul * np.array(f4p5["/suscept_func/IM_DENSITY"])
rechi_m4p5 = vert_mul * np.array(f4p5["/suscept_func/RE_MAGNETIC"])
imchi_m4p5 = vert_mul * np.array(f4p5["/suscept_func/IM_MAGNETIC"])

fdim_bos4p5       = bosgrid4p5.shape[0]
mom_dim_bos4p5    = rechi_t4p5.shape[1]
ffactor_count4p5  = rechi_t4p5.shape[2]

#BETA=5

bosgrid5 = np.array(f5["/suscept_func/bgrid"])
rechi_t5 = vert_mul * np.array(f5["/suscept_func/RE_TRIPLET"])
imchi_t5 = vert_mul * np.array(f5["/suscept_func/IM_TRIPLET"])
rechi_s5 = vert_mul * np.array(f5["/suscept_func/RE_SINGLET"])
imchi_s5 = vert_mul * np.array(f5["/suscept_func/IM_SINGLET"])
rechi_d5 = vert_mul * np.array(f5["/suscept_func/RE_DENSITY"])
imchi_d5 = vert_mul * np.array(f5["/suscept_func/IM_DENSITY"])
rechi_m5 = vert_mul * np.array(f5["/suscept_func/RE_MAGNETIC"])
imchi_m5 = vert_mul * np.array(f5["/suscept_func/IM_MAGNETIC"])

fdim_bos5       = bosgrid5.shape[0]
mom_dim_bos5    = rechi_t5.shape[1]
ffactor_count5  = rechi_t5.shape[2]
#==================================================================================================
#
#       PLOT N1: SUSCEPTIBILITIES IN THE DIFFERENT CHANNELS AS A FUNCTIONS OF TEMPERATURES
#
#==================================================================================================


omega = 0

omega1 = omega + (fdim_bos1 - 1)/2 
omega2 = omega + (fdim_bos2 - 1)/2 
omega4 = omega + (fdim_bos4 - 1)/2 

print omega2, omega4

#q = (0,0)

if (PATCH_COUNT1 == 16):
    q1 = 5
elif (PATCH_COUNT1 == 4):
    q1 = 0
elif (PATCH_COUNT1== 64):
    q1 = 27

if (PATCH_COUNT2 == 16):
    q2 = 5
elif (PATCH_COUNT2 == 4):
    q2 = 0
elif (PATCH_COUNT2== 64):
    q2 = 27

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#SINGLET
rechis_s = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#DENSITY
rechid_s = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#MAGNETIC
rechim_s = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create s*-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_ss = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_ss = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#SINGLET                                                                                                    
rechis_ss = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_ss = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#DENSITY                                1                                            1                       
rechid_ss = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_ss = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC   

rechim_ss = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_ss = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create d-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_d = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_d = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#SINGLET 

rechis_d = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_d = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#DENSITY  

rechid_d = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_d = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC                              1                                            1                        
rechim_d = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_d = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create px-wave array for different temperatures
FF1 = 3
FF2 = 3 

#TRIPLET
rechit_px = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_px = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#SINGLET                                                       
rechis_px = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_px = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#DENSITY                                                       
rechid_px = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_px = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#MAGNETIC                                                      
rechim_px = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_px = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create py-wave array for different temperatures
FF1 = 4
FF2 = 4 

#TRIPLET
rechit_py = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_py = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                
#SINGLET                                                       
rechis_py = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_py = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#DENSITY                                                       
rechid_py = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_py = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#MAGNETIC                                                      
rechim_py = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_py = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


beta_array = np.array([1./BETA5,1./BETA4p5,1./BETA4,1./BETA2, 1./BETA1])


#--- Helper functions

def plotchi( use_pl, arrs, arrss, arrd, arrpx, arrpy, legend):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrss,"ro", label = "$s^{*}$")
    pl.plot( beta_array, arrd, "gs", label = "$d_{x^2-y^2}$")
    pl.plot( beta_array, arrpx,"m^", label = "$p_{x}$")
    pl.plot( beta_array, arrpy,"kv", label = "$p_{y}$")
    pl.grid(linestyle='dotted')
    pl.xlim([0.0, 1.5])
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    
    return




fig = pl.figure(figsize=cm2inch(15.0,15.0))

#pl.subplots_adjust( wspace=0.6, hspace=15)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega_n=$" + str(omega)+ r"             $q = (0,0)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_ss, rechit_d, rechit_px, rechit_py, False) 
pl.ylabel(RE + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_ss, rechis_d, rechis_px, rechis_py, True) 
pl.ylabel(RE + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_ss, rechid_d, rechid_px, rechid_py, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_ss, rechim_d, rechim_px, rechim_py, False) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()
#pl.tight_layout()

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# RE PLOT FOR POSTER! Density, Magnetic, SC = 1/2(s+t)

#pl.subplots_adjust( wspace=0.4, hspace=0.4)

fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_ss+rechis_ss),0.5*(rechit_d+rechis_d), 0.5*(rechit_px+rechis_px), 0.5*(rechit_py+ rechis_py), False) 
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.ylabel(RE + r"\chi^{SC}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_ss, rechid_d, rechid_px, rechid_py, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_ss, rechim_d, rechim_px, rechim_py, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()


pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

#pl.subplots_adjust( wspace=0.4, hspace=0.4)
fig = pl.figure(figsize=cm2inch(15.0,15.0))
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_ss, imchit_d, imchit_px, imchit_py, False) 
pl.ylabel(IM + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_ss, imchis_d, imchis_px, imchis_py, True) 
pl.ylabel(IM + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_ss, imchid_d, imchid_px, imchid_py, False) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_ss, imchim_d, imchim_px, imchim_py, False) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#====================================================================================
#
#                                   q = (pi,pi)
#
#====================================================================================

if (PATCH_COUNT1 == 16):
    q1 = 15
elif (PATCH_COUNT1 == 4):
    q1 = 3
elif (PATCH_COUNT1== 64):
    q1 = 63

if (PATCH_COUNT2 == 16):
    q2 = 15
elif (PATCH_COUNT2 == 4):
    q2 = 3
elif (PATCH_COUNT2== 64):
    q2 = 63

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                              
#SINGLET                                                       
rechis_s = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#DENSITY                                                      
rechid_s = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                               
#MAGNETIC                                                     
rechim_s = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create s*-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_ss = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_ss = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#SINGLET                                1                                                                    
rechis_ss = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_ss = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#DENSITY                                1                                                                   
rechid_ss = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_ss = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC   

rechim_ss = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_ss = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create d-wave array for different temperatures
FF1 = 1
FF2 = 2 

#TRIPLET
rechit_d = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_d = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#SINGLET 

rechis_d = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_d = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#DENSITY  

rechid_d = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_d = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#MAGNETIC                              1                                            1                        
rechim_d = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_d = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
#Create p-wave(symmetric) array for different temperatures
FF1 = 3
FF2 = 4 

#TRIPLET
rechit_ps = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_ps = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#SINGLET                                                                                                    
rechis_ps = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_ps = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#DENSITY                                                                                                   
rechid_ps = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_ps = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#MAGNETIC                                                                                                  
rechim_ps = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_ps = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] + imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] + imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create p-wave(antisymmetric) array for different temperatures
FF1 = 3
FF2 = 4 

#TRIPLET
rechit_pa = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_pa = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#SINGLET                                                                                                    
rechis_pa = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_pa = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#DENSITY                                                                                                    
rechid_pa = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_pa = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                            
#MAGNETIC                                                                                                   
rechim_pa = np.array([rechi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_pa = np.array([imchi_m5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF1, 0, 0, 0, 0] - imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create coskx-wave array for different temperatures
FF1 = 1
FF2 = 1 

#TRIPLET
rechit_cx = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_cx = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                
#SINGLET                                                        
rechis_cx = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_cx = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                        
#DENSITY                                                                                                
rechid_cx = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_cx = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                
#MAGNETIC                                                       
rechim_cx = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_cx = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


#Create cosky-wave array for different temperatures
FF1 = 2
FF2 = 2 

#TRIPLET
rechit_cy = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_cy = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#SINGLET                                                                                                  
rechis_cy = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_cy = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#DENSITY                                                                                                  
rechid_cy = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_cy = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#MAGNETIC                                                                                                 
rechim_cy = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_cy = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])


beta_array = np.array([1./BETA5, 1./BETA4p5,1./BETA4,1./BETA2, 1./BETA1])

#--- Helper functions

def plotchi( use_pl, arrs, arrcx, arrcy, arrps, arrpa, legend):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrcx,"ro", label = "$s*$")
    pl.plot( beta_array, arrcy,"gs", label = "$d_{x^2-y^2}$")
    pl.plot( beta_array, arrps,"m^", label = "$p_{x}+p_y$")
    pl.plot( beta_array, arrpa,"kv", label = "$p_{x}-p_y$")
    pl.grid(linestyle='dotted')
    pl.xlim([0.0, 1.5])
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    
    return


fig = pl.figure(figsize=cm2inch(15.0,15.0))

#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_ss, rechit_d, rechit_ps, rechit_pa, False) 
pl.ylabel(RE + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_ss, rechis_d, rechis_ps, rechis_pa, True) 
pl.ylabel(RE + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_ss, rechid_d, rechid_ps, rechid_pa, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_ss, rechim_d, rechim_ps, rechim_pa, False) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# RE PLOT POSTER!
fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_cx+rechis_cx), 0.5*(rechit_cy+ rechis_cy), 0.5*(rechit_ps+rechis_ps), 0.5*(rechit_pa+rechis_pa), False) 
pl.ylabel(RE + r"\chi^{SC}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_ss, rechid_d, rechid_ps, rechid_pa, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_ss, rechim_d, rechim_ps, rechim_pa, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#IM

fig = pl.figure(figsize=cm2inch(15.0,15.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (\pi,\pi)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_ss, imchit_d, imchit_ps, imchit_pa, False) 
pl.ylabel(IM + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_ss, imchis_d, imchis_ps, imchis_pa, True) 
pl.ylabel(IM + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_ss, imchid_d, imchid_ps, imchid_pa, False) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_ss, imchim_d, imchim_ps, imchim_pa, False) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#===================================================================================================================================
#
#                                       q = (0,pi)
#
#==================================================================================================================================

#q=(0,pi)

if (PATCH_COUNT1 == 16):
    q1 = 7
elif (PATCH_COUNT1 == 4):
    q1 = 1
elif (PATCH_COUNT1== 64):
    q1 = 31

if (PATCH_COUNT2 == 16):
    q2 = 7
elif (PATCH_COUNT2 == 4):
    q2 = 1
elif (PATCH_COUNT2== 64):
    q2 = 31

#Create s-wave array for different temperatures
FF1 = 0
FF2 = 0 

#TRIPLET
rechit_s = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_s = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                         
#SINGLET                                                                                                 
rechis_s = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_s = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                         
#DENSITY                                                                                                 
rechid_s = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_s = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                         
#MAGNETIC                                                                                                
rechim_s = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_s = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create cos(kx)-wave array for different temperatures

FF1 = 1
FF2 = 1 

#TRIPLET
rechit_cx = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_cx = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#SINGLET                                                                                                  
rechis_cx = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_cx = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                         
#DENSITY                                                                                                  
rechid_cx = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_cx = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#MAGNETIC                                                                                                 
rechim_cx = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_cx = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create cos(ky)-wave array for different temperatures

FF1 = 2
FF2 = 2 

#TRIPLET
rechit_cy = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_cy = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#SINGLETy                                                                                                 
rechis_cy = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_cy = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#DENSITY                                                                                                  
rechid_cy = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_cy = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#MAGNETIY                                                                                                 
rechim_cy = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_cy = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create sin(kx)-wave array for different temperatures

FF1 = 3
FF2 = 3 

#TRIPLET
rechit_sx = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_sx = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#SINGLET                                                                                                  
rechis_sx = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_sx = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#DENSITY                                                                                                  
rechid_sx = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_sx = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#MAGNETIC                                                                                                 
rechim_sx = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_sx = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])

#Create sin(ky)-wave array for different temperatures

FF1 = 4
FF2 = 4 

#TRIPLET
rechit_sy = np.array([rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchit_sy = np.array([imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_t2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_t1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#SINGLETy                                                                                                 
rechis_sy = np.array([rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchis_sy = np.array([imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_s2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_s1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#DENSITY                                                                                                  
rechid_sy = np.array([rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchid_sy = np.array([imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_d1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
                                                                                                          
#MAGNETIY                                                                                                 
rechim_sy = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], rechi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])
imchim_sy = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4p5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q2, FF1, FF2, 0, 0, 0, 0], imchi_m1[omega1, q2, FF1, FF2, 0, 0, 0, 0]])



beta_array = np.array([1./BETA5, 1./BETA4p5,1./BETA4,1./BETA2, 1./BETA1])

#--- Helper functions

def plotchi( use_pl, arrs, arrcx, arrcy, arrsx, arrsy, legend):
    pl.plot( beta_array, arrs, "bx", label = "$s$")
    pl.plot( beta_array, arrcx,"ro", label = "$\cos(k_x)$")
    pl.plot( beta_array, arrcy,"gs", label = "$\cos(k_y)}$")
    pl.plot( beta_array, arrsx,"m^", label = "$\sin(k_x)$")
    pl.plot( beta_array, arrsy,"kv", label = "$\sin(k_y)$")
    pl.grid(linestyle='dotted')
    pl.xlim([0.0, 1.5])
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    
    return


fig = pl.figure(figsize=cm2inch(15.0,15.0))

#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(2,2,1), rechit_s, rechit_cx, rechit_cy, rechit_sx, rechit_sy, False) 
pl.ylabel(RE + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), rechis_s, rechis_cx, rechis_cy, rechis_sx, rechis_sy, True) 
pl.ylabel(RE + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), rechid_s, rechid_cx, rechid_cy, rechid_sx, rechid_sy, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), rechim_s, rechim_cx, rechim_cy, rechim_sx, rechim_sy, False) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_re_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

# PLOT RE POSTER!
fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(1,3,1), 0.5*(rechit_s+rechis_s), 0.5*(rechit_cx+rechis_cx), 0.5*(rechit_cy+rechis_cy), 0.5*(rechit_sx+rechis_sx), 0.5*(rechit_sy+rechis_sy), False) 
pl.ylabel(RE + r"\chi^{t}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,2), rechid_s, rechid_cx, rechid_cy, rechid_sx, rechid_sy, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(1,3,3), rechim_s, rechim_cx, rechim_cy, rechim_sx, rechim_sy, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_poster_re_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#IM

fig = pl.figure(figsize=cm2inch(15.0,15.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,\pi)$" )

plotchi( pl.subplot(2,2,1), imchit_s, imchit_cx, imchit_cy, imchit_sx, imchit_sy, False) 
pl.ylabel(IM + r"\chi^{t}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,2), imchis_s, imchis_cx, imchis_cy, imchis_sx, imchis_sy, True) 
pl.ylabel(IM + r"\chi^{s}$")
#pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,3), imchid_s, imchid_cx, imchid_cy, imchid_sx, imchid_sy, False) 
pl.ylabel(IM + r"\chi^{d}$")
pl.xlabel(r"$T/t$") 
plotchi( pl.subplot(2,2,4), imchim_s, imchim_cx, imchim_cy, imchim_sx, imchim_sy, False) 
pl.ylabel(IM + r"\chi^{m}$")
pl.xlabel(r"$T/t$") 
pl.tight_layout()

pl.savefig("plots/suscept_im_Om="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()



#===============================================================================================================
#
#                                       SELF-ENERGY
#
#===============================================================================================================


fname2 = "dat/carsten_comparison/U2/PREPROC/dat_U2_TP0_Beta2_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_maxkpos4_allsymm.h5"  
fname4 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta4_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  
fname4p5 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta4_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  
fname5 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta5_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  


if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname4 = str(sys.argv[1])

if len(sys.argv) > 1:
    fname4p5 = str(sys.argv[1])

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname4= fname4.rstrip('\n') # strip newline of fname
f4 = h5py.File(fname4, "r")

fname4p5 = fname4p5.rstrip('\n') # strip newline of fname
f4p5 = h5py.File(fname4p5, "r")

fname5 = fname5.rstrip('\n') # strip newline of fname
f5 = h5py.File(fname5, "r")

#--- Read
#BETA=2
siggrid2  = f2["/Sig/fgrid"]
resig2    = prefact * np.array(f2["/Sig/RE"])
imsig2    = prefact * np.array(f2["/Sig/IM"])
mom_dim2  = resig2.shape[1]
fdim2     = resig2.shape[0]
fdimo22   = fdim2/2

#BETA=4
siggrid4  = f4["/Sig/fgrid"]
resig4    = prefact * np.array(f4["/Sig/RE"])
imsig4    = prefact * np.array(f4["/Sig/IM"])
mom_dim4  = resig4.shape[1]
fdim4     = resig4.shape[0]
fdimo24   = fdim4/2

#BETA=4.5
siggrid4p5  = f4p5["/Sig/fgrid"]
resig4p5    = prefact * np.array(f4p5["/Sig/RE"])
imsig4p5    = prefact * np.array(f4p5["/Sig/IM"])
mom_dim4p5  = resig4p5.shape[1]
fdim4p5     = resig4p5.shape[0]
fdimo24p5   = fdim4p5/2

#BETA=5
siggrid5  = f5["/Sig/fgrid"]
resig5    = prefact * np.array(f5["/Sig/RE"])
imsig5    = prefact * np.array(f5["/Sig/IM"])
mom_dim5  = resig5.shape[1]
fdim5     = resig5.shape[0]
fdimo25   = fdim5/2

#==============================================================
#
#       PLOT: SE_REAL AND SE_IMAG FOR DIFFERENT K-POINTS
#
#==============================================================

#--- Helper functions
 #q(0,0)
if PATCH_COUNT1==4:
    q1= 0
elif PATCH_COUNT1==16:
    q1 = 5
elif PATCH_COUNT1==64:
    q1 = 27

if PATCH_COUNT2==4:
    q2 = 0
elif PATCH_COUNT2==16:
    q2 = 5
elif PATCH_COUNT2==64:
    q2 = 27

def plotSig( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx',  label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go',  label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv',  label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms',  label=r'$\beta=5$')
    pl.xlim([min(siggrid2),max(siggrid2)])
    pl.grid(linestyle='dotted')
    use_pl.set_title(string)
    pl.legend(loc=1, prop={'size':6})
    return

def plotSigre( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx', label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go', label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv', label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms',  label=r'$\beta=5$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    return

fig = pl.figure(figsize=cm2inch(14.0,8.0))

pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (0,0)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig2[:,q2,0,0], resig4[:,q1,0,0], resig4p5[:,q1,0,0], resig5[:,q1,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig2[:,q2,0,0], imsig4[:,q1,0,0], imsig4p5[:,q1,0,0], imsig5[:,q1,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--- Helper functions
 #q=(0,pi)

if PATCH_COUNT1==4:
    q1 = 1
elif PATCH_COUNT1==16:
    q1 = 7
elif PATCH_COUNT1==64:
    q1 = 31

if PATCH_COUNT2==4:
    q2 = 1
elif PATCH_COUNT2==16:
    q2 = 7
elif PATCH_COUNT2==64:
    q2 = 31

def plotSig( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx',label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go',label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv', label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.xlim([min(siggrid2),max(siggrid2)])
    pl.grid(linestyle='dotted')
    use_pl.set_title(string)
    pl.legend(loc=1, prop={'size':6})
    return

def plotSigre( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx',  label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go',  label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv', label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    return

fig = pl.figure(figsize=cm2inch(14.0,8.0))
pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (0,\pi)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig2[:,q2,0,0], resig4[:,q1,0,0], resig4p5[:,q1,0,0], resig5[:,q1,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig2[:,q2,0,0], imsig4[:,q1,0,0], imsig4p5[:,q1,0,0], imsig5[:,q1,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Helper functions
 #q=(pi/2,pi/2)

if PATCH_COUNT1==4:
    print "NO (pi/2,pi/2) point for PATCH_COUNT=4"
elif PATCH_COUNT1==16:
    q1 = 10
elif PATCH_COUNT1==64:
    q1 = 45

if PATCH_COUNT2==4:
    print "NO (pi/2,pi/2) point for PATCH_COUNT=4"
elif PATCH_COUNT2==16:
    q2 = 10
elif PATCH_COUNT2==64:
    q2 = 45


def plotSig( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx', label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go', label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv',  label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.xlim([min(siggrid2),max(siggrid2)])
    pl.grid(linestyle='dotted')
    use_pl.set_title(string)
    pl.legend(loc=1, prop={'size':6})
    return

def plotSigre( use_pl, arr1, arr2, arr3,arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx', label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go', label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv', label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    return

fig = pl.figure(figsize=cm2inch(14.0,8.0))
pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (\pi/2,\pi/2)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig2[:,q2,0,0], resig4[:,q1,0,0], resig4p5[:,q1,0,0], resig5[:,q1,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig2[:,q2,0,0], imsig4[:,q1,0,0], imsig4p5[:,q1,0,0], imsig5[:,q1,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()


#--- Save to file
pl.savefig("plots/Sig_k=(pi2,pi2).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--- Helper functions
 #q=(pi,pi)

if PATCH_COUNT1==4:
    q1 = 3
elif PATCH_COUNT1==16:
    q1 = 15
elif PATCH_COUNT1==64:
    q1 = 63

if PATCH_COUNT2==4:
    q2 = 3
elif PATCH_COUNT2==16:
    q2 = 15
elif PATCH_COUNT2==64:
    q2 = 63

def plotSig( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx', label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go', label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv',  label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    pl.legend(loc=1, prop={'size':6})
    return

def plotSigre( use_pl, arr1, arr2, arr3, arr4, string ):
    pl.plot( siggrid2, arr1[:], 'bx', label=r"$\beta=2$")
    pl.plot( siggrid4, arr2[:], 'go', label=r'$\beta=4$')
    pl.plot( siggrid4p5, arr3[:], 'rv', label=r'$\beta=4.5$')
    pl.plot( siggrid5, arr4[:], 'ms', label=r'$\beta=5$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid2),max(siggrid2)])
    use_pl.set_title(string)
    return

fig = pl.figure(figsize=cm2inch(14.0,8.0))
pl.suptitle(r"$U=$" + str(UINT1) + r"     $q = (\pi,\pi)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig2[:,q2,0,0], resig4[:,q1,0,0], resig4p5[:,q1,0,0], resig5[:,q1,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig2[:,q2,0,0], imsig4[:,q1,0,0], imsig4p5[:,q1,0,0], imsig5[:,q1,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#=========================================================================================================
#
#                              FITTING OF THE LINEAR REGIME OF THE SELF-ENERGY
#                                      -> QP WEIGHT Z TO BE CALCULATED
#
#                           WARNING: TOO HIGH TEMPERATURE..FITTING JUST TWO POINTS!!!
#========================================================================================================

freqset4   = np.array([(2*i+1)*pi/4. for i in range(0,2)])
freqset4p5 = np.array([(2*i+1)*pi/4.5 for i in range(0,2)])
freqset5 = np.array([(2*i+1)*pi/5. for i in range(0,2)])

#q(0,0)
if PATCH_COUNT1==4:
    q = 0
elif PATCH_COUNT1==16:
    q = 5
elif PATCH_COUNT1==64:
    q = 27

#alpha4_00 = (imsig4[0+fdimo24,q,0,0]-imsig4[1+fdimo24,q,0,0])/(freqset4[0]-freqset4[1])
#alpha4p5_00 = (imsig4p5[0+fdimo24p5,q,0,0]-imsig4p5[1+fdimo24p5,q,0,0])/(freqset4p5[0]-freqset4p5[1])
#alpha5_00 = (imsig5[0+fdimo25,q,0,0]-imsig5[1+fdimo25,q,0,0])/(freqset5[0]-freqset5[1])

alpha4_00 = (-imsig4[0+fdimo24,q,0,0])/(pi/4.0)
alpha4p5_00 = (-imsig4p5[0+fdimo24p5,q,0,0])/(pi/4.5)
alpha5_00 = (-imsig5[0+fdimo25,q,0,0])/(pi/5.0)

Z4_00   = 1.0/(alpha4_00+1)
Z4p5_00 = 1.0/(alpha4p5_00+1)
Z5_00 = 1.0/(alpha5_00+1)


 #q=(0,pi)

if PATCH_COUNT1==4:
    q = 1
elif PATCH_COUNT1==16:
    q = 7
elif PATCH_COUNT1==64:
    q = 31

alpha4_0pi = (-imsig4[0+fdimo24,q,0,0])/(pi/4.0)
alpha4p5_0pi = (-imsig4p5[0+fdimo24p5,q,0,0])/(pi/4.5)
alpha5_0pi = (-imsig5[0+fdimo25,q,0,0])/(pi/5.0)
#alpha4_0pi   = (imsig4[0+fdimo24,q,0,0]-imsig4[1+fdimo24,q,0,0])/(freqset4[0]-freqset4[1])
#alpha4p5_0pi = (imsig4p5[0+fdimo24p5,q,0,0]-imsig4p5[1+fdimo24p5,q,0,0])/(freqset4p5[0]-freqset4p5[1])
#alpha5_0pi = (imsig5[0+fdimo25,q,0,0]-imsig5[1+fdimo25,q,0,0])/(freqset5[0]-freqset5[1])

Z4_0pi   = 1.0/(alpha4_0pi+1)
Z4p5_0pi = 1.0/(alpha4p5_0pi+1)
Z5_0pi = 1.0/(alpha5_0pi+1)

print Z5_0pi

#q=(pi/2,pi/2)

if PATCH_COUNT1==4:
    print "NO (pi/2,pi/2) point for PATCH_COUNT=4"
elif PATCH_COUNT1==16:
    q = 10
elif PATCH_COUNT1==64:
    q = 45

alpha4_pi2pi2 = (-imsig4[0+fdimo24,q,0,0])/(pi/4.0)
alpha4p5_pi2pi2 = (-imsig4p5[0+fdimo24p5,q,0,0])/(pi/4.5)
alpha5_pi2pi2 = (-imsig5[0+fdimo25,q,0,0])/(pi/5.0)
#alpha4_pi2pi2   = (imsig4[0+fdimo24,q,0,0]-imsig4[1+fdimo24,q,0,0])/(freqset4[0]-freqset4[1])
#alpha4p5_pi2pi2 = (imsig4p5[0+fdimo24p5,q,0,0]-imsig4p5[1+fdimo24p5,q,0,0])/(freqset4p5[0]-freqset4p5[1])
#alpha5_pi2pi2 = (imsig5[0+fdimo25,q,0,0]-imsig5[1+fdimo25,q,0,0])/(freqset5[0]-freqset5[1])

Z4_pi2pi2   = 1.0/(alpha4_pi2pi2+1)
Z4p5_pi2pi2 = 1.0/(alpha4p5_pi2pi2+1)
Z5_pi2pi2 = 1.0/(alpha5_pi2pi2+1)
 
#q=(pi,pi)

if PATCH_COUNT1==4:
    q = 3
elif PATCH_COUNT1==16:
    q = 15
elif PATCH_COUNT1==64:
    q = 63


alpha4_pipi = (-imsig4[0+fdimo24,q,0,0])/(pi/4.0)
alpha4p5_pipi = (-imsig4p5[0+fdimo24p5,q,0,0])/(pi/4.5)
alpha5_pipi = (-imsig5[0+fdimo25,q,0,0])/(pi/5.0)
#alpha4_pipi   = (imsig4[0+fdimo24,q,0,0]-imsig4[1+fdimo24,q,0,0])/(freqset4[0]-freqset4[1])
#alpha4p5_pipi = (imsig4p5[0+fdimo24p5,q,0,0]-imsig4p5[1+fdimo24p5,q,0,0])/(freqset4p5[0]-freqset4p5[1])
#alpha5_pipi   = (imsig5[0+fdimo25,q,0,0]-imsig5[1+fdimo25,q,0,0])/(freqset5[0]-freqset5[1])

print alpha4_pipi, alpha4p5_pipi

Z4_pipi   = 1.0/(alpha4_pipi+1)
Z4p5_pipi = 1.0/(alpha4p5_pipi+1)
Z5_pipi = 1.0/(alpha5_pipi+1)

x = [0,1,2,3]
x_plot = [-1,0,1,2,3,4]
labels=["", "$(0,0)$", "$(0,\pi)$", "$(\pi/2,\pi/2)$", "$(\pi,\pi)$", ""]

def plotZ( use_pl):
    pl.plot( x, np.array([Z4_00, Z4_0pi, Z4_pi2pi2, Z4_pipi]), 'go', label="$beta=4$")
    #pl.plot( x, np.array([Z4p5_00, Z4p5_0pi, Z4p5_pi2pi2, Z4p5_pipi]), 'rv', label="$beta=4.5$")
    pl.plot( x, np.array([Z5_00, Z5_0pi, Z5_pi2pi2, Z5_pipi]), 'ms', label="$beta=5$")
    pl.xticks(x_plot, labels)
    #pl.ylim([1.0, 1.0009])
    pl.ylabel("Z")
    pl.grid(linestyle='dotted')
    pl.legend(loc=4, prop={'size':8})
    return

#pl.suptitle(r"$U=$" + str(UINT1) )

fig = pl.figure(figsize=cm2inch(15.0,15.0))
#--- Plot physical
plotZ( pl.subplot(1,1,1) ) 

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Z.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

