#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: cx219
"""
from wobster.WOBSTER import *
from wobster.w90io import *

#---------------------------------------------
# input parameters:
num_bands = 8       #number of wannier functions
num_kpoints = 18**3 #total number of K-points
energy_min = -15    #DOS energy minimum [eV]
energy_max = 20     #DOS energy maximum [eV]
NEDOS = 200         #Total number of points between energy_min and energy_max
SIGMA = 0.2         #smearing factor

# read-in data
eigenvals = read_eig("wannier90.eig",8,num_kpoints)
U_matrix,kpoints = read_u_matrix("wannier90_u.mat",num_bands,num_bands,num_kpoints)
hop,Rlatt,deg = read_hr()

#---------------------------------------------
# get me DOS
# **Note that spin degeneracy is not considered here**
print "Calculating total dos."
dos = get_dos(eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('DOS_total.dat', dos)
#----plot total dos----
fig, ax = plt.subplots()
ax.plot(dos[:,0],dos[:,1])
fig.savefig("total_DOS.png", dpi=300)
print "done."

#---------------------------------------------
# get me projected DOS (WOOP)
# we want all orbitals to be within origin cell
R1 = [0,0,0]
R2 = [0,0,0]

# For s-s between Si atoms
#---------------------------------------------
print "Calculating WOOPs for s orbital."
dos_s = get_WOOP(U_matrix,kpoints, R1,R2,0,0, eigenvals,num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('DOS_s.dat', dos_s)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_s[:,0],dos_s[:,1])
fig.savefig("WOOP_s.png", dpi=300)
print "done."


# p
print "Calculating WOOPs for p orbitals."
dos_px = get_WOOP(U_matrix, kpoints, R1,R2,5,5, eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('DOS_px.dat', dos_px)
dos_py = get_WOOP(U_matrix, kpoints, R1,R2,6,6, eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('DOS_py.dat', dos_py)
dos_pz = get_WOOP(U_matrix, kpoints, R1,R2,7,7, eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('DOS_pz.dat', dos_pz)
dos_p = dos_px
dos_p[:,1] = dos_px[:,1]+dos_py[:,1]+dos_pz[:,1]
np.savetxt('DOS_p.dat', dos_p)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_p[:,0],dos_p[:,1])
fig.savefig("WOOP_p.png", dpi=300)
print "done."
#---------------------------------------------

#---------------------------------------------
# WOHP calculations
#---------------------------------------------

# For s-s between Si atoms
#---------------------------------------------
print "Calculating WOHPs between s-s orbitals."

## we are calculating wohp between R1_wan1 and R2_wan2
wan1 = 0
wan2 = 4
# the meaningful number for R is the difference between R1 and R2
latt_diff = np.array(R2, dtype=int) - np.array(R1, dtype=int)
# now calculate the index for cell R
cell_indx_hop = np.where([np.all(latt==latt_diff) for latt in Rlatt])[0][0]
# using the index, get me the hopping.
hopping = hop[cell_indx_hop,wan1,wan2].real

dos_ss = get_WOHP(hopping,U_matrix,kpoints,R1,R2,wan1,wan2,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_ss.dat', dos_ss)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_ss[:,0],dos_ss[:,1])
fig.savefig("WOHP_ss.png", dpi=300)
print "done."

# For s-p between Si atoms
#---------------------------------------------
print "Calculating WOHPs between s-p orbitals."
# the meaningful number for R is the difference between R1 and R2
latt_diff = np.array(R2, dtype=int) - np.array(R1, dtype=int)
# now calculate the index for cell R
cell_indx_hop = np.where([np.all(latt==latt_diff) for latt in Rlatt])[0][0]

# s-px
dos_spx = get_WOHP(hop[cell_indx_hop,0,5].real,U_matrix,kpoints,R1,R2,0,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_spx.dat', dos_spx)
# s-py
dos_spy = get_WOHP(hop[cell_indx_hop,0,6].real,U_matrix,kpoints,R1,R2,0,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_spy.dat', dos_spy)
# s-pz
dos_spz = get_WOHP(hop[cell_indx_hop,0,7].real,U_matrix,kpoints,R1,R2,0,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_spz.dat', dos_spz)
dos_sp = dos_spx.copy()
dos_sp[:,1] = dos_spx[:,1]+dos_spy[:,1]+dos_spz[:,1]
np.savetxt('WOHP_sp.dat', dos_sp)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_sp[:,0],dos_sp[:,1])
fig.savefig("WOHP_sp.png", dpi=300)
print "done."

# For p-p between Si atoms
#---------------------------------------------
print "Calculating WOHPs between p-p orbitals."
# the meaningful number for R is the difference between R1 and R2
latt_diff = np.array(R2, dtype=int) - np.array(R1, dtype=int)
# now calculate the index for cell R
cell_indx_hop = np.where([np.all(latt==latt_diff) for latt in Rlatt])[0][0]

dos_pxpx = get_WOHP(hop[cell_indx_hop,1,5].real,U_matrix,kpoints,R1,R2,1,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpx.dat', dos_pxpx)
dos_pxpy = get_WOHP(hop[cell_indx_hop,1,6].real,U_matrix,kpoints,R1,R2,1,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpy.dat', dos_pxpy)
dos_pxpz = get_WOHP(hop[cell_indx_hop,1,7].real,U_matrix,kpoints,R1,R2,1,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pxpz.dat', dos_pxpz)
dos_pypx = get_WOHP(hop[cell_indx_hop,2,5].real,U_matrix,kpoints,R1,R2,2,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypx.dat', dos_pypx)
dos_pypy = get_WOHP(hop[cell_indx_hop,2,6].real,U_matrix,kpoints,R1,R2,2,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypy.dat', dos_pypy)
dos_pypz = get_WOHP(hop[cell_indx_hop,2,7].real,U_matrix,kpoints,R1,R2,2,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pypz.dat', dos_pypz)
dos_pzpx = get_WOHP(hop[cell_indx_hop,3,5].real,U_matrix,kpoints,R1,R2,3,5,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpx.dat', dos_pzpx)
dos_pzpy = get_WOHP(hop[cell_indx_hop,3,6].real,U_matrix,kpoints,R1,R2,3,6,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpy.dat', dos_pzpy)
dos_pzpz = get_WOHP(hop[cell_indx_hop,3,7].real,U_matrix,kpoints,R1,R2,3,7,eigenvals, num_kpoints, energy_min, energy_max, NEDOS, SIGMA)
np.savetxt('WOHP_pzpz.dat', dos_pzpz)
dos_pp = dos_pxpx
dos_pp[:,1] = dos_pxpx[:,1]+dos_pxpy[:,1]+dos_pxpz[:,1]+dos_pypx[:,1]+dos_pypy[:,1]+dos_pypz[:,1]+dos_pzpx[:,1]+dos_pzpy[:,1]+dos_pzpz[:,1]
np.savetxt('WOHP_pp.dat', dos_pp)
#----plot partial dos----
fig, ax = plt.subplots()
ax.plot(dos_pp[:,0],dos_pp[:,1])
fig.savefig("WOHP_pp.png", dpi=300)
print "done."

#---------------------------------------------
