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

fname1 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta1_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  
fname2 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta2_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  

fname4 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta4_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  
fname5 = "dat/carsten_comparison/U2/PREPROC/OMFL/dat_U2_Beta5_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_ALLSYMM_FNN_mkp2.h5"  

if len(sys.argv) > 1:
    fname4 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname5 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")


fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")

fname4 = fname4.rstrip('\n') # strip newline of fname
f4 = h5py.File(fname4, "r")


fname5 = fname5.rstrip('\n') # strip newline of fname
f5 = h5py.File(fname5, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals1 = f1["/Params"].attrs.values()
parVals2 = f2["/Params"].attrs.values()
parVals4 = f4["/Params"].attrs.values()
parVals5 = f5["/Params"].attrs.values()

## BETA = 1
UINT1 =  parVals1[0] # follows order in output.cpp
BETA1 =  parVals1[1]
B1    =  parVals1[2]
MU1   =  parVals1[3]

# BETA = 2
UINT2 =  parVals2[0] # follows order in output.cpp
BETA2 =  parVals2[1]
B2    =  parVals2[2]
MU2   =  parVals2[3]

# BETA = 4
UINT4 =  parVals4[0] # follows order in output.cpp
BETA4 =  parVals4[1]
B4 =     parVals4[2]
MU4 =     parVals4[3]

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

vert_mul = 1.0/4.0/pi/pi
#vert_mul = 1.0
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

rechi_sc1 = vert_mul * np.array(f1["/susc_func/RE_SC"])
imchi_sc1 = vert_mul * np.array(f1["/susc_func/IM_SC"])
rechi_d1  = vert_mul * np.array(f1["/susc_func/RE_D"])
imchi_d1  = vert_mul * np.array(f1["/susc_func/IM_D"])
rechi_m1  = vert_mul * np.array(f1["/susc_func/RE_M"])
imchi_m1  = vert_mul * np.array(f1["/susc_func/IM_M"])

fdim_bos1       = rechi_sc1.shape[0]
mom_dim_bos1    = rechi_sc1.shape[1]
ffactor_count1  = rechi_sc1.shape[2]

#BETA=2

rechi_sc2 = vert_mul * np.array(f2["/susc_func/RE_SC"])
imchi_sc2 = vert_mul * np.array(f2["/susc_func/IM_SC"])
rechi_d2  = vert_mul * np.array(f2["/susc_func/RE_D"])
imchi_d2  = vert_mul * np.array(f2["/susc_func/IM_D"])
rechi_m2  = vert_mul * np.array(f2["/susc_func/RE_M"])
imchi_m2  = vert_mul * np.array(f2["/susc_func/IM_M"])

fdim_bos2       = rechi_sc2.shape[0]
mom_dim_bos2    = rechi_sc2.shape[1]
ffactor_count2  = rechi_sc2.shape[2]

#BETA=4

rechi_sc4 = vert_mul * np.array(f4["/susc_func/RE_SC"])
imchi_sc4 = vert_mul * np.array(f4["/susc_func/IM_SC"])
rechi_d4  = vert_mul * np.array(f4["/susc_func/RE_D"])
imchi_d4  = vert_mul * np.array(f4["/susc_func/IM_D"])
rechi_m4  = vert_mul * np.array(f4["/susc_func/RE_M"])
imchi_m4  = vert_mul * np.array(f4["/susc_func/IM_M"])

fdim_bos4       = rechi_sc4.shape[0]
mom_dim_bos4    = rechi_sc4.shape[1]
ffactor_count4  = rechi_sc4.shape[2]

#BETA=5

rechi_sc5 = vert_mul * np.array(f5["/susc_func/RE_SC"])
imchi_sc5 = vert_mul * np.array(f5["/susc_func/IM_SC"])
rechi_d5  = vert_mul * np.array(f5["/susc_func/RE_D"])
imchi_d5  = vert_mul * np.array(f5["/susc_func/IM_D"])
rechi_m5  = vert_mul * np.array(f5["/susc_func/RE_M"])
imchi_m5  = vert_mul * np.array(f5["/susc_func/IM_M"])

fdim_bos5       = rechi_sc5.shape[0]
mom_dim_bos5    = rechi_sc5.shape[1]
ffactor_count5  = rechi_sc5.shape[2]
#
##BETA=4
#
#bosgrid4 = np.array(f4["/suscept_func/bgrid"])
#rechi_t4 = vert_mul * np.array(f4["/suscept_func/RE_TRIPLET"])
#imchi_t4 = vert_mul * np.array(f4["/suscept_func/IM_TRIPLET"])
#rechi_s4 = vert_mul * np.array(f4["/suscept_func/RE_SINGLET"])
#imchi_s4 = vert_mul * np.array(f4["/suscept_func/IM_SINGLET"])
#rechi_d4 = vert_mul * np.array(f4["/suscept_func/RE_DENSITY"])
#imchi_d4 = vert_mul * np.array(f4["/suscept_func/IM_DENSITY"])
#rechi_m4 = vert_mul * np.array(f4["/suscept_func/RE_MAGNETIC"])
#imchi_m4 = vert_mul * np.array(f4["/suscept_func/IM_MAGNETIC"])
#
#fdim_bos4       = bosgrid4.shape[0]
#mom_dim_bos4    = rechi_t4.shape[1]
#ffactor_count4  = rechi_t4.shape[2]
#
##BETA=5
#
#bosgrid5 = np.array(f5["/suscept_func/bgrid"])
#rechi_t5 = vert_mul * np.array(f5["/suscept_func/RE_TRIPLET"])
#imchi_t5 = vert_mul * np.array(f5["/suscept_func/IM_TRIPLET"])
#rechi_s5 = vert_mul * np.array(f5["/suscept_func/RE_SINGLET"])
#imchi_s5 = vert_mul * np.array(f5["/suscept_func/IM_SINGLET"])
#rechi_d5 = vert_mul * np.array(f5["/suscept_func/RE_DENSITY"])
#imchi_d5 = vert_mul * np.array(f5["/suscept_func/IM_DENSITY"])
#rechi_m5 = vert_mul * np.array(f5["/suscept_func/RE_MAGNETIC"])
#imchi_m5 = vert_mul * np.array(f5["/suscept_func/IM_MAGNETIC"])
#
#fdim_bos5       = bosgrid5.shape[0]
#mom_dim_bos5    = rechi_t5.shape[1]
#ffactor_count5  = rechi_t5.shape[2]
#
#==================================================================================================
#
#       PLOT N1: SUSCEPTIBILITIES IN THE DIFFERENT CHANNELS AS A FUNCTIONS OF TEMPERATURES
#
#==================================================================================================


omega = 0

omega1 = omega + (fdim_bos1 - 1)/2 
omega2 = omega + (fdim_bos2 - 1)/2 
omega4 = omega + (fdim_bos4 - 1)/2 
omega4 = omega + (fdim_bos5 - 1)/2 

print omega4

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

#Create d-wave SC array for different temperatures
FF1 = 1
FF2 = 2 

#SC
rechisc_d = np.array([rechi_sc5[omega4, q1, FF1, FF1, 0, 0, 0, 0]-rechi_sc5[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_sc4[omega4, q1, FF1, FF1, 0, 0, 0, 0]-rechi_sc4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_sc2[omega2, q1, FF1, FF1, 0, 0, 0, 0]-rechi_sc2[omega2, q1, FF1, FF2, 0, 0, 0, 0], rechi_sc1[omega1, q1, FF1, FF1, 0, 0, 0, 0]-rechi_sc1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ] )
imchisc_d = np.array([imchi_sc5[omega4, q1, FF1, FF1, 0, 0, 0, 0], imchi_sc4[omega4, q1, FF1, FF1, 0, 0, 0, 0], imchi_sc2[omega2, q1, FF1, FF1, 0, 0, 0, 0], imchi_sc1[omega1, q1, FF1, FF1, 0, 0, 0, 0] ] )

##TRIPLET
#rechit_d = np.array([rechi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0] ])
#imchit_d = np.array([imchi_t5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_t4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_t4[omega4, q1, FF1, FF2, 0, 0, 0, 0] ])
#
##SINGLET 
#rechis_d = np.array([rechi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0] ])
#imchis_d = np.array([imchi_s5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_s4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_s4[omega4, q1, FF1, FF2, 0, 0, 0, 0] ])

rechisc_carst_d = np.array([ 0.660558,  0.570955, 0.396861, 0.325328, 0.241404])

#Create d-wave Pomeranchuk array for different temperatures
FF1 = 1
FF2 = 2 

#DENSITY  
rechipom_d = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],rechi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])
imchipom_d = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],imchi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])

rechipom_carst_d = np.array([ 0.436934,  0.397436, 0.305684, 0.261817, 0.205879])

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

#Create AFM array for different temperatures
FF1 = 0
FF2 = 0 
                                                
rechim_s = np.array([rechi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_m2[omega2, q1, FF1, FF2, 0, 0, 0, 0],rechi_m1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])
imchim_s = np.array([imchi_m5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_m4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_m2[omega2, q1, FF1, FF2, 0, 0, 0, 0],imchi_m1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])

rechim_carst_s = np.array([ 1.80095, 1.29725, 0.67271, 0.499891, 0.332516])

#Create d-wave DENSITY array for different temperatures
FF1 = 1
FF2 = 2 

rechid_d = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],rechi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])
imchid_d = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],imchi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])

#Create p-wave DENSITY array for different temperatures
FF1 = 3
FF2 = 4 

rechid_p = np.array([rechi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],rechi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], rechi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],rechi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - rechi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])
imchid_p = np.array([imchi_d5[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d5[omega4, q1, FF1, FF2, 0, 0, 0, 0],imchi_d4[omega4, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d4[omega4, q1, FF1, FF2, 0, 0, 0, 0], imchi_d2[omega2, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d2[omega2, q1, FF1, FF2, 0, 0, 0, 0],imchi_d1[omega1, q1, FF1, FF1, 0, 0, 0, 0] - imchi_d1[omega1, q1, FF1, FF2, 0, 0, 0, 0] ])

rechid_carst_d = np.array([ 0.660558,  0.570955, 0.396861, 0.325328, 0.241404])

beta_array = np.array([1./BETA5,1./BETA4, 1./BETA2, 1./BETA1])
beta_array_carst = np.array([1./BETA5, 0.255, 0.44, 0.58, 0.86])


#--- Helper functions

def plotchi( use_pl, arrSCd, arrSCdcar, arrPomd, arrPomdcar, arrAF, arrAFcar, arrDd, arrDdcar,arrDp, legend):
    pl.plot( beta_array, arrSCd, "ro-", label = r"$\chi_{dSC}^{-1}$")
    pl.plot( beta_array_carst, arrSCdcar,"ro--", mfc='none',label = r"$\chi_{dSC, c}^{-1}$")
    pl.plot( beta_array, arrPomd, "gs-", label = r"$\chi_{dPom}^{-1}$")
    pl.plot( beta_array_carst, arrPomdcar,"gs--",mfc='none', label = r"$\chi_{dPom,c}^{-1}$")
    pl.plot( beta_array, arrAF,"bd-", label = r"$\chi_{AF}^{-1}$")
    pl.plot( beta_array_carst, arrAFcar,"bd--", mfc='none',label = r"$\chi_{AF,c}^{-1}$")
    pl.plot( beta_array, arrDd,"m^-", label = "$\chi_{dDW}^{-1}$")
    pl.plot( beta_array_carst, arrDdcar,"m^--", mfc='none', label = r"$\chi_{dDW,c}^{-1}$")
    pl.plot( beta_array, arrDp,"kd-", label = "$\chi_{pDW}^{-1}$")
    pl.grid(linestyle='dotted')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 5.0])
    if(legend):
        pl.legend(loc=1, prop={'size':8})
    
    return




# RE PLOT FOR POSTER! Density, Magnetic, SC = 1/2(s+t)

#pl.subplots_adjust( wspace=0.4, hspace=0.4)

fig = pl.figure(figsize=cm2inch(15.0,15.0))
#pl.suptitle(r"$U=$" + str(UINT1) + r"     $\Omega=$" + str(omega)+ r"     $q = (0,0)$" )
plotchi( pl.subplot(1,1,1), 1.0/(rechisc_d), 1.0/(rechisc_carst_d), 1.0/(rechipom_d), 1.0/(rechipom_carst_d), 1.0/(rechim_s), 1.0/(rechim_carst_s), 1.0/(rechid_d), 1.0/(rechid_carst_d),1.0/rechid_p, True) 
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.xlabel(r"$T/t$") 
pl.tight_layout()
#ax.grid()

pl.savefig("plots/suscept_comparison_carsten_data.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

