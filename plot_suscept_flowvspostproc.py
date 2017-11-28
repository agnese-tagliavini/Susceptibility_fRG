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

#fname1 = "dat/carsten_comparison/U2/suscept/suscept_U2_TP0_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  
#fname2 = "dat/carsten_comparison/U2/suscept_flow/dat_U2_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_mkp4_allsymm_withSE_with_4pi2.h5"  
#fname2 = "dat/dat_U2_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_mkp2_allsymm_with_4pi2_check_triasympt_new.h5"  
#fname1 = "dat/suscept_U2_Beta5_PFCB256_HUB_SU2_2D_INTFL.h5"  

fname2 = "dat/RES_FUNC_SCHEMES/SCHEME_4/dat_U2_Beta5_PFCB128_HUB_OMFL_KAT_2LOOPRESFUNC_SU2_2D_4PISYMM_mkp2_ALLSYMM_WITHSE.h5"  
fname1 = "dat/RES_FUNC_SCHEMES/suscept/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_OMFL_2LOOP.h5"  

#fname2 = "dat/check_bare_bubble_INTFL/dat_U2_Beta5_PFCB128_HUB_OMFL_SU2_2D_4PISYMM_mkp2_allsymm_noSE_with_4pi2_swave_bare_withasympt_novert.h5"  
#fname1 = "dat/check_bare_bubble_INTFL/suscept_U2_Beta5_PFCB128_HUB_SU2_2D_bare_withasympt.h5" 

#fname2 = "dat/check_INTFL/dat_U2_Beta5_PFCB128_HUB_INTFL_SU2_2D_4PISYMM_ALLSYMM_LOCAL_mkp2.h5"  
#fname1 = "dat/check_INTFL/suscept_U2_TP0_Beta5_PFCB128_HUB_SU2_2D_4PISYMM_mkp4_allsymm_withSE.h5"  

if len(sys.argv) > 1:
    fname1 = str(sys.argv[1])
if len(sys.argv) > 1:
    fname2 = str(sys.argv[1])

fname1 = fname1.rstrip('\n') # strip newline of fname
f1 = h5py.File(fname1, "r")

fname2 = fname2.rstrip('\n') # strip newline of fname
f2 = h5py.File(fname2, "r")


os.system('mkdir -p log')
os.system('mkdir -p plots')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals1 = f1["/Params"].attrs.values()
parVals2 = f2["/Params"].attrs.values()

# BETA = 5 post-processing

UINT1 =  parVals1[0] # follows order in output.cpp
BETA1 =  parVals1[1]
B1 =     parVals1[2]
MU1 =     parVals1[3]

# BETA = 5 during the flow

UINT2 =  parVals2[0] # follows order in output.cpp
BETA2 =  parVals2[1]
B2 =     parVals2[2]
MU2 =     parVals2[3]

T_PRIME= 0.0

MAX_KPOS = 2
PATCH_COUNT=4*MAX_KPOS*MAX_KPOS

pi = math.pi

prefact = 1.0

vert_mul = 1.0/4.0/pi/pi
#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
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

bosgrid_plot = np.array(f1["/suscept_func/bgrid"])
rechi_t1 =  np.array(f1["/suscept_func/RE_TRIPLET"])
imchi_t1 =  np.array(f1["/suscept_func/IM_TRIPLET"])
rechi_s1 =  np.array(f1["/suscept_func/RE_SINGLET"])
imchi_s1 =  np.array(f1["/suscept_func/IM_SINGLET"])
rechi_d1 =  np.array(f1["/suscept_func/RE_DENSITY"])
imchi_d1 =  np.array(f1["/suscept_func/IM_DENSITY"])
rechi_m1 =  np.array(f1["/suscept_func/RE_MAGNETIC"])
imchi_m1 =  np.array(f1["/suscept_func/IM_MAGNETIC"])

fdim_bos1       = bosgrid_plot.shape[0]
mom_dim_bos1    = rechi_t1.shape[1]
ffactor_count1  = rechi_t1.shape[2]

#BETA=5 tprime = 0.3

#bosgrid_plot = np.array(f2["/suscept_func/bgrid"])
rechi_sc2 =  vert_mul*np.array(f2["/susc_func/RE_SC"])
imchi_sc2 =  vert_mul*np.array(f2["/susc_func/IM_SC"])
rechi_d2 =   vert_mul*np.array(f2["/susc_func/RE_D"])
imchi_d2 =   vert_mul*np.array(f2["/susc_func/IM_D"])
rechi_m2 =   vert_mul*np.array(f2["/susc_func/RE_M"])
imchi_m2 =   vert_mul*np.array(f2["/susc_func/IM_M"])

fdim_bos2       = rechi_sc2.shape[0]
mom_dim_bos2    = rechi_sc2.shape[1]
ffactor_count2  = rechi_sc2.shape[2]

fdim_arr = np.array([fdim_bos1,fdim_bos2])

fdim_min = fdim_arr.min()

bosgrid_plot = np.array([(2*n)*pi/BETA1 for n in range(-(fdim_min-1)/2, (fdim_min-1)/2)])

#==================================================================================================
#
#       PLOT N1: SUSCEPTIBILITIES IN THE DIFFERENT CHANNELS AS A FUNCTIONS OF TEMPERATURES
#
#==================================================================================================


omega = 0

omega1 = omega + (fdim_bos1 - 1)/2 
omega2 = omega + (fdim_bos2 - 1)/2 


#q = (0,0)

if (PATCH_COUNT == 16):
    q = 5
elif (PATCH_COUNT == 4):
    q = 0
elif (PATCH_COUNT== 64):
    q = 27

FF1 = 0
FF2 = 0

def plotchifreq( use_pl, arr1, arr2, legend):
    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $s$")
    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $s$")
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    return


fig = pl.figure(figsize=cm2inch(21.0,10.0))
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,0)$      $\beta =$"+ str(BETA1))

plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
pl.ylabel(RE + r"\chi^{sc}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_swave_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#FF1 = 1
#FF2 = 2
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF1, 0,0,0,0]-arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $d_{x^2-y^2}$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF1, 0,0,0,0]-arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $d_{x^2-y^2}$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,0)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_dwave_q=(0,0).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF1, 0,0,0,0]+arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $s^{*}$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF1, 0,0,0,0]+arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $s^{*}$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,0)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_sswave_q=(0,0).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()

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

FF1 = 0
FF2 = 0

def plotchifreq( use_pl, arr1, arr2,legend):
    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $s$")
    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $s$")
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    return


fig = pl.figure(figsize=cm2inch(21.0,10.0))
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (\pi,\pi)$      $\beta =$"+ str(BETA1))

plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
pl.ylabel(RE + r"\chi^{sc}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_swave_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#FF1 = 1
#FF2 = 2
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF1, 0,0,0,0]-arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $d_{x^2-y^2}$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF1, 0,0,0,0]-arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $d_{x^2-y^2}$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (\pi,\pi)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_dwave_q=(pi,pi).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF1, 0,0,0,0]+arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $s^{*}$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF1, 0,0,0,0]+arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $s^{*}$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (\pi,\pi)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_sswave_q=(pi,pi).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#

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


FF1 = 0
FF2 = 0

def plotchifreq( use_pl, arr1, arr2, legend):
    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $s$")
    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF1, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $s$")
    if(legend):
        pl.legend(loc=1, prop={'size':6})
    return


fig = pl.figure(figsize=cm2inch(21.0,10.0))
pl.subplots_adjust( wspace=0.4, hspace=0.4)
pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,\pi)$      $\beta =$"+ str(BETA1))

plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
pl.ylabel(RE + r"\chi^{sc}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
pl.ylabel(RE + r"\chi^{d}$")
pl.xlabel(r"$\Omega$") 
plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
pl.ylabel(RE + r"\chi^{m}$")
pl.xlabel(r"$\Omega$") 

pl.savefig("plots/suscept_re_swave_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#FF1 = 1
#FF2 = 2
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF1, FF1, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $\cos(k_y)$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2, q, FF1, FF1, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $\cos(k_x)$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,\pi)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_coskx_q=(0,pi).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()
#
#def plotchifreq( use_pl, arr1, arr2, legend):
#    pl.plot( bosgrid_plot, np.array([ arr1[ n+(fdim_bos1-1)/2, q, FF2, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "bx", label = "POSTPROC, $\cos(k_y)$")
#    pl.plot( bosgrid_plot, np.array([ arr2[ n+(fdim_bos2-1)/2,q, FF2, FF2, 0,0,0,0] for n in range(-(fdim_min-1)/2,(fdim_min-1)/2)]), "ro", label = "FLOW, $\cos(k_y)$")
#    if(legend):
#        pl.legend(loc=1, prop={'size':6})
#    return
#
#
#fig = pl.figure(figsize=cm2inch(21.0,10.0))
#pl.subplots_adjust( wspace=0.4, hspace=0.4)
#pl.suptitle(r"$U=$" + str(UINT1) +r" $q = (0,\pi)$      $\beta =$"+ str(BETA1))
#
#plotchifreq( pl.subplot(1,3,1), 0.5*(rechi_t1+rechi_s1), rechi_sc2, False) 
#pl.ylabel(RE + r"\chi^{sc}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,2), rechi_d1, rechi_d2, False) 
#pl.ylabel(RE + r"\chi^{d}$")
#pl.xlabel(r"$\Omega$") 
#plotchifreq( pl.subplot(1,3,3), rechi_m1, rechi_m2, True) 
#pl.ylabel(RE + r"\chi^{m}$")
#pl.xlabel(r"$\Omega$") 
#
#pl.savefig("plots/suscept_re_cosky_q=(0,pi).png", dpi=150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()




