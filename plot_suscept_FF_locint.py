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


fname = "dat/carsten_comparison/U2/suscept/suscept_U2_Beta4p5_PFCB16_HUB_SU2_2D_4PISYMM_maxkpos4_allsymm_withSE.h5"  
#fname = "dat/" + most_recently_edited

if len(sys.argv) > 1:
    fname = str(sys.argv[1])

fname = fname.rstrip('\n') # strip newline of fname
f = h5py.File(fname, "r")

os.system('mkdir -p log')
os.system('mkdir -p plots')

#os.system('rm plots/prev_Vert.png 2> log/plot.log')
#os.system('rm plots/prev_phi.png 2> log/plot.log')
#os.system('rm plots/prev_P_func.png 2> log/plot.log')
#os.system('rm plots/prev_R_func.png 2> log/plot.log')

#os.system('mv plots/Vert.png plots/prev_Vert.png 2> log/plot.log')
#os.system('mv plots/phi.png plots/prev_phi.png 2> log/plot.log')
#os.system('mv plots/P_func.png plots/prev_P_func.png 2> log/plot.log')
#os.system('mv plots/R_func.png plots/prev_R_func.png 2> log/plot.log')


#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

parVals = f["/Params"].attrs.values()

UINT =  parVals[0] # follows order in output.cpp
BETA =  parVals[1]
B =     parVals[2]
MU =     parVals[3]
#T_PRIME=parVals[4]

T_PRIME=  0 #parVals[4]

MAX_KPOS = 4
PATCH_COUNT=4*MAX_KPOS*MAX_KPOS
#PATCH_COUNT=16   # if patch_count=16, add dif and inv for Simple_patching is used (predefined matrices)
shift=0

nu=0
nup=0

FF1=0
FF2=1


pi = math.pi
#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
#pl.rc('text', usetex=True)
#pl.rc('text.latex', preamble='\usepackage{amsmath}')

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
    form_factor_arr = np.array([1//math.sqrt(4*pi*pi), 
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
bosgrid = np.array(f["/suscept_func/bgrid"])
rechi_t = vert_mul * np.array(f["/suscept_func/RE_TRIPLET"])
imchi_t = vert_mul * np.array(f["/suscept_func/IM_TRIPLET"])
rechi_s = vert_mul * np.array(f["/suscept_func/RE_SINGLET"])
imchi_s = vert_mul * np.array(f["/suscept_func/IM_SINGLET"])
rechi_d = vert_mul * np.array(f["/suscept_func/RE_DENSITY"])
imchi_d = vert_mul * np.array(f["/suscept_func/IM_DENSITY"])
rechi_m = vert_mul * np.array(f["/suscept_func/RE_MAGNETIC"])
imchi_m = vert_mul * np.array(f["/suscept_func/IM_MAGNETIC"])

fdim_bos = bosgrid.shape[0]
mom_dim_bos = rechi_t.shape[1]
ffactor_count = rechi_t.shape[2]

freq_bos_grid = np.array([2*i*pi/BETA for i in range(-(fdim_bos-1)/4, (fdim_bos-1)/4+1)])

y_range_fact = 1.0

#--- Helper functions

def plotchi( use_pl, arr, string ):
#    use_pl.set_aspect(1.0)
    pl.pcolormesh( np.arange(mom_dim_bos+1),freq_bos_grid, np.array([arr[i+(fdim_bos-1)/2,:,FF1,FF2,0,0,0,0] for i in range(-(fdim_bos-1)/4, (fdim_bos-1)/4+1)]))
    pl.ylim([y_range_fact*min(freq_bos_grid),y_range_fact*max(freq_bos_grid)])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return



pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}$" )
pl.ylabel(r"$\Omega$")
pl.xlabel(r"$Q$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}$" )
pl.xlabel(r"$Q$")
pl.tight_layout()

pl.savefig("plots/Suscept_re_FF1="+str(FF1)+"_FF2="+str(FF2)+".png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), imchi_t, IM + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), imchi_s, IM + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), imchi_d, IM + r"\chi^{d}$" )
pl.ylabel(r"$\Omega$")
pl.xlabel(r"$Q$")
plotchi( pl.subplot(2,2,4), imchi_m, IM + r"\chi^{m}$" )
pl.xlabel(r"$Q$")

pl.tight_layout()

pl.savefig("plots/Suscept_im_FF1="+str(FF1)+"_FF2="+str(FF2)+".png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#---------------------------------Chi LOCAL-------------------------------------------
#--- Helper functions

def plotchi( use_pl, arr, string ):
#    use_pl.set_aspect(1.0)
    pl.plot( freq_bos_grid, np.array([np.sum(arr[i+(fdim_bos-1)/2,:,0,0,0,0,0,0], axis=0)/(mom_dim_bos) for i in range(-(fdim_bos-1)/4, (fdim_bos-1)/4+1)]))
    pl.xlim([y_range_fact*min(bosgrid),y_range_fact*max(bosgrid)])
    use_pl.set_title( string , fontsize=10 )
    return

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}_{loc}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}_{loc}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}_{loc}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}_{loc}$" )
pl.xlabel(r"$\Omega$")

pl.tight_layout()

pl.savefig("plots/Suscept_re_loc.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), imchi_t, IM + r"\chi^{t}_{loc}$" )
plotchi( pl.subplot(2,2,2), imchi_s, IM + r"\chi^{s}_{loc}$" )
plotchi( pl.subplot(2,2,3), imchi_d, IM + r"\chi^{d}_{loc}$" )
pl.xlabel(r"$\Omega$")
plotchi( pl.subplot(2,2,4), imchi_m, IM + r"\chi^{m}_{loc}$" )
pl.xlabel(r"$\Omega$")

pl.tight_layout()

pl.savefig("plots/Suscept_im_loc.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


#--------------------------------------Chi PLOTTING MATRIX FOR FIXED OMEGA AND Q ------------------------------------------

print("Plotting suceptibilities ...")

omega = 0

#q = (0,0)

if (PATCH_COUNT == 16):
    q = 5
elif (PATCH_COUNT == 4):
    q = 0

elif (PATCH_COUNT== 64):
    q = 27


#--- Helper functions

def plotchi( use_pl, arr, string ):
#    use_pl.set_aspect(1.0)
    pl.pcolormesh( np.arange(ffactor_count+1),np.arange(ffactor_count+1), arr[omega+(fdim_bos-1)/2,q,:,:,0,0,0,0])
    use_pl.set_title( string , fontsize=10 )
    pl.colorbar(shrink=0.6) 
    return



pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")
pl.tight_layout()

pl.savefig("plots/Suscept_re_W="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), imchi_t, IM + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), imchi_s, IM + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), imchi_d, IM + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), imchi_m, IM + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")

pl.tight_layout()

pl.savefig("plots/Suscept_im_W="+str(omega)+"_q=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#q=(pi,pi)
if (PATCH_COUNT == 16):
    q = 15
elif (PATCH_COUNT == 4):
    q = 3

elif (PATCH_COUNT== 64):
    q = 63

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")
pl.tight_layout()

pl.savefig("plots/Suscept_re_W="+str(omega)+"_q=(pi,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#q=(pi/2,pi/2)
if (PATCH_COUNT == 16):
    q = 15
elif (PATCH_COUNT == 4):
    q = 3

elif (PATCH_COUNT== 64):
    q = 45

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")
pl.tight_layout()

pl.savefig("plots/Suscept_re_W="+str(omega)+"_q=(pi2,pi2).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
#IM

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), imchi_t, IM + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), imchi_s, IM + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), imchi_d, IM + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), imchi_m, IM + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")

pl.tight_layout()

pl.savefig("plots/Suscept_im_W="+str(omega)+"_q=(pi2,pi2).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#q=(0,pi)
if (PATCH_COUNT == 16):
    q = 7
elif (PATCH_COUNT == 4):
    q = 1

elif (PATCH_COUNT== 64):
    q = 31

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), rechi_t, RE + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), rechi_s, RE + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), rechi_d, RE + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), rechi_m, RE + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")
pl.tight_layout()

pl.savefig("plots/Suscept_re_W="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#IM

pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) )

plotchi( pl.subplot(2,2,1), imchi_t, IM + r"\chi^{t}$" )
plotchi( pl.subplot(2,2,2), imchi_s, IM + r"\chi^{s}$" )
plotchi( pl.subplot(2,2,3), imchi_d, IM + r"\chi^{d}$" )
pl.ylabel(r"$f_1$")
pl.xlabel(r"$f_2$")
plotchi( pl.subplot(2,2,4), imchi_m, IM + r"\chi^{m}$" )
pl.xlabel(r"$f_2$")

pl.tight_layout()

pl.savefig("plots/Suscept_im_W="+str(omega)+"_q=(0,pi).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()
