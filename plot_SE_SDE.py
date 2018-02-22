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


fname = "/home/agnese/Coding/fRG_code/2D/fRGdyn/dat/dat_U2_Beta2_PFCB128_MKP2_HUB_OMFL_KAT_2LOOP_MLOOP_NUMLOOP1_SECORR_SU2_2D_4PISYMM_ALLSYMM_WITHSE.h5"  
fnameaf = "dat/suscept_U2_Beta2_PFCB128_HUB_KAT_2LOOP_SU2_2D_3LOOP_OMFL.h5"  


f = h5py.File(fname, "r")
faf = h5py.File(fnameaf, "r")


parVals = f["/Params"].attrs.values()

UINT =  parVals[0] # follows order in output.cpp
BETA =  parVals[1]
B =     parVals[2]
MU =    parVals[3]

MAX_KPOS = 2
PATCH_COUNT =4*MAX_KPOS*MAX_KPOS

pi = math.pi

vert_mul = 1.0/4.0/pi/pi
prefact = 1.0
#prefact1 = 4.0 * pi * pi
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

siggrid  = f["/Sig/fgrid"]
resig    = prefact * np.array(f["/Sig/RE"])
imsig    = prefact * np.array(f["/Sig/IM"])
mom_dim  = resig.shape[1]
fdim     = resig.shape[0]
fdimo2   = fdim/2

siggridaf  = faf["/Sig/fgrid"]
resigaf    = prefact * np.array(faf["/Sig/RE"])
imsigaf    = prefact * np.array(faf["/Sig/IM"])
mom_dimaf  = resigaf.shape[1]
fdimaf     = resigaf.shape[0]
fdimo2af   = fdimaf/2
 
 #q(0,0)
if PATCH_COUNT==4:
    q= 0
elif PATCH_COUNT==16:
    q = 5
elif PATCH_COUNT==64:
    q = 27

def plotSig( use_pl, arr1, arr2, string ):
    pl.plot( siggrid, arr1[:], 'bx',  label=r"$\Sigma$")
    pl.plot( siggridaf, arr2[:], 'go',  label=r'$\Sigma_{SD}$')
    pl.xlim([min(siggrid),max(siggrid)])
    pl.grid(linestyle='dotted')
    use_pl.set_title(string)
    pl.legend(loc=1, prop={'size':6})
    return

def plotSigre( use_pl, arr1, arr2,string ):
    pl.plot( siggrid, arr1[:], 'bx', label=r"$\Sigma$")
    pl.plot( siggridaf, arr2[:], 'go', label=r'$\Sigma_{SD}$')
    pl.grid(linestyle='dotted')
    pl.xlim([min(siggrid),max(siggrid)])
    use_pl.set_title(string)
    return

fig = pl.figure(figsize=cm2inch(14.0,8.0))

pl.suptitle(r"$U=$" + str(UINT) + r"     $q = (0,0)$" )

#--- Plot physical
plotSigre( pl.subplot(1,2,1), resig[:,q,0,0], resigaf[:,q,0,0], RE + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
plotSig( pl.subplot(1,2,2), imsig[:,q,0,0], imsigaf[:,q,0,0], IM + "\Sigma(i\omega)$" ) 
pl.xlabel("$\omega$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig_k=(0,0).png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()


