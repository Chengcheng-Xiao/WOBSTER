#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Python class to get TB hamiltonian from realspace wannier functions.
@author: cx219
"""
from wobster.WOBSTER import *

#---------------------------------------
# some naming convention
# U_matrix[k_num,num_bands,num_wann,dat]
# U_matrix_ct[k,num_wann,num_bands,dat]
# eigenvals[kpoint_number,band_number]
# kpoints[kpoint_number,dir]

# input parameters:
num_bands = 8
num_kpoints = 18**3
energy_min = -15
energy_max = 20
NEDOS = 200
SIGMA = 0.2

# generate stuff
eigenvals = get_eig("wannier90.eig",8,num_kpoints)
U_matrix,kpoints = get_u_matrix("wannier90_u.mat",num_bands,num_bands,num_kpoints)

# # get me DOS
# dos = get_dos(eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('DOS_total.dat', dos)
# #----plot total dos----
# fig, ax = plt.subplots()
# ax.plot(dos[:,0],dos[:,1])
# fig.savefig("total_DOS.png", dpi=300)


#----plot total dos----
# fig, ax = plt.subplots()
# ax.plot(dos[:,0],dos[:,1])
# fig.savefig("total_dos.png", dpi=300)


#---------------------------------------------
# get me projected DOS (WOOP)
# we want all orbitals to be within origin cell
R1 = [0,0,0]
R2 = [0,0,0]

# # s-s
# dos_s = get_WOOP(U_matrix,kpoints,R1,0,eigenvals,num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('DOS_s.dat', dos_s)
# #----plot partial dos----
# fig, ax = plt.subplots()
# ax.plot(dos_s[:,0],dos_s[:,1])
# fig.savefig("WOOP_s.png", dpi=300)


# # p
# dos_px = get_WOOP(U_matrix,kpoints,R1,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('DOS_px.dat', dos_px)
# dos_py = get_WOOP(U_matrix,kpoints,R1,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('DOS_py.dat', dos_py)
# dos_pz = get_WOOP(U_matrix,kpoints,R1,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('DOS_pz.dat', dos_pz)
# dos_p = dos_px
# dos_p[:,1] = dos_px[:,1]+dos_py[:,1]+dos_pz[:,1]
# np.savetxt('DOS_p.dat', dos_p)
# #----plot partial dos----
# fig, ax = plt.subplots()
# ax.plot(dos_p[:,0],dos_p[:,1])
# fig.savefig("WOOP_p.png", dpi=300)
#---------------------------------------------

#---------------------------------------------
# get me WOHP

# # s-s
# dos_ss = get_WOHP(-3.055005,U_matrix,kpoints,R1,R2,0,4,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('WOHP_ss.dat', dos_ss)
# #----plot partial dos----
# fig, ax = plt.subplots()
# ax.plot(dos_ss[:,0],dos_ss[:,1])
# fig.savefig("WOHP_ss.png", dpi=300)

# # s-p
# dos_spx = get_WOHP(-1.701416,U_matrix,kpoints,R1,R2,0,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('WOHP_spx.dat', dos_spx)
# dos_spy = get_WOHP(-4.251791,U_matrix,kpoints,R1,R2,0,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('WOHP_spy.dat', dos_spy)
# dos_spz = get_WOHP(-2.415476,U_matrix,kpoints,R1,R2,0,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
# np.savetxt('WOHP_spz.dat', dos_spz)
# dos_sp = dos_spx
# dos_sp[:,1] = dos_spx[:,1]+dos_spy[:,1]+dos_spz[:,1]
# np.savetxt('WOHP_sp.dat', dos_sp)
# #----plot partial dos----
# fig, ax = plt.subplots()
# ax.plot(dos_sp[:,0],dos_sp[:,1])
# fig.savefig("WOHP_sp.png", dpi=300)

# p-p
dos_pxpx = get_WOHP(-1.348359,U_matrix,kpoints,R1,R2,1,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpx.dat', dos_pxpx)
dos_pxpy = get_WOHP( 2.028108,U_matrix,kpoints,R1,R2,1,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpy.dat', dos_pxpy)
dos_pxpz = get_WOHP( 1.197320,U_matrix,kpoints,R1,R2,1,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpz.dat', dos_pxpz)
dos_pypx = get_WOHP( 2.028108,U_matrix,kpoints,R1,R2,2,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypx.dat', dos_pypx)
dos_pypy = get_WOHP( 2.714971,U_matrix,kpoints,R1,R2,2,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypy.dat', dos_pypy)
dos_pypz = get_WOHP( 2.867865,U_matrix,kpoints,R1,R2,2,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypz.dat', dos_pypz)
dos_pzpx = get_WOHP( 1.197320,U_matrix,kpoints,R1,R2,3,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpx.dat', dos_pzpx)
dos_pzpy = get_WOHP( 2.867865,U_matrix,kpoints,R1,R2,3,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpy.dat', dos_pzpy)
dos_pzpz = get_WOHP(-0.511690,U_matrix,kpoints,R1,R2,3,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpz.dat', dos_pzpz)
dos_pp = dos_pxpx
dos_pp[:,1] = dos_pxpx[:,1]+dos_pxpy[:,1]+dos_pxpz[:,1]+dos_pypx[:,1]+dos_pypy[:,1]+dos_pypz[:,1]+dos_pzpx[:,1]+dos_pzpy[:,1]+dos_pzpz[:,1]
np.savetxt('WOHP_pp.dat', dos_pp)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_pp[:,0],dos_pp[:,1])
fig.savefig("WOHP_pp.png", dpi=300)

#---------------------------------------------
