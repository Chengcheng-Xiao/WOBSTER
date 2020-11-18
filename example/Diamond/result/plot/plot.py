import numpy as np
import matplotlib as mpl
# mpl.use('agg')
mpl.rcParams['figure.dpi']= 300
import matplotlib.pyplot as plt

# font stuff
font = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 12}

font1 = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 20}

font2 = {'family' : 'DejaVu Sans',
        'weight' : 'light',
        'size'   : 8}

mpl.rc('font', **font)

COHP = np.loadtxt('COHPCAR.lobster')
WOHP_sp = np.loadtxt('WOHP_sp.dat')
WOHP_ss = np.loadtxt('WOHP_ss.dat')

# account for spin degeneracy
WOHP_ss[:,1] *= 2
WOHP_sp[:,1] *= 2


E_fermi = 9.2
fig, axs = plt.subplots(1,2,figsize=(10,3))
axs[0].plot(WOHP_ss[:,0]-E_fermi,WOHP_ss[:,1],label='s-s')
axs[0].plot(WOHP_sp[:,0]-E_fermi,WOHP_sp[:,1],label='s-p')
axs[0].hlines(0,WOHP_ss[:,0].min()-E_fermi,WOHP_ss[:,0].max()-E_fermi,linewidth=0.5, color='k',linestyles='dashed')
axs[0].vlines(0,(-COHP_sp).min(),(-COHP_sp).max(),linewidth=0.5, color='k',linestyles='dashed')
axs[0].set_xlim(WOHP_ss[:,0].min()-E_fermi,WOHP_ss[:,0].max()-E_fermi)
axs[0].set_ylim((-COHP_sp).min(),(-COHP_sp).max())
axs[0].set_xlabel("Energy (eV)",**font1)
axs[0].set_title("WOHP",**font1)
axs[0].legend()

COHP_sp = COHP[:,7]+COHP[:,9]+COHP[:,11]
axs[1].plot(COHP[:,0],-COHP[:,5],label='s-s')
axs[1].plot(COHP[:,0],-(COHP_sp),label='s-p')
axs[1].hlines(0,COHP[:,0].min(),COHP[:,0].max(),linewidth=0.5, color='k',linestyles='dashed')
axs[1].vlines(0,(-COHP_sp).min(),(-COHP_sp).max(),linewidth=0.5, color='k',linestyles='dashed')
axs[1].set_xlim(COHP[:,0].min(),COHP[:,0].max())
axs[1].set_ylim((-COHP_sp).min(),(-COHP_sp).max())
axs[1].set_xlabel("Energy (eV)",**font1)
axs[1].set_title("COHP",**font1)
axs[1].legend()
