import numpy as np
import matplotlib.pyplot as plt

def coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt):
    '''
    Calculate projection coefficient:

        coeff = exp(i*2pi*K*R)*\tilde{U^K_{ibnd, iwan}}

    input:
        U_matrix: full unitary matrix
        ibnd: band number
        iwan: wannier function number
        R: lattice vectore, should be array(3)
        kpoints: array contains all kpoints info
        ikpt: number kpoint

    return: coef
    Format: coef
    '''
    coef = np.conj(complex(U_matrix[ikpt,ibnd,iwan,0],U_matrix[ikpt,ibnd,iwan,1]))
    kdR = 2*np.pi*np.dot(kpoints[ikpt],R)
    #coef *= ((2*np.pi)**3/V)*np.exp(complex(0,kdR))
    coef *= np.exp(complex(0,kdR))
    return coef


def get_dos(eigenvals, num_kpoints, e_min, e_max, nedos, T=0.1):
    '''
    calculate density of states by:

    DOS(\epsilon) = \sum_j \delta(E_j-\epsilon)

    \delta function is expressed as:
    \delta = \frac{1}{\sqrt{\pi*T}}*exp(-(E_j-\epsilon)**2/T)

    where T is the semaring factor.

    ** Note that spin degeneracy is not considered here. **

    returns: dos
    format:  dos[ieng,energy/dos]
    '''
    dos = np.zeros((nedos,2), dtype=float)
    energy = np.linspace(e_min,e_max,nedos)
    for ieng in range(nedos):
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                dos[ieng,0] = energy[ieng]
                dos[ieng,1] += 1/np.sqrt(np.pi*T)*np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    return dos

# def get_proj_dos(coeff,iwan,eigenvals, num_kpoints, e_min, e_max, nedos, T=0.1):
#     dos = np.zeros((nedos,2), dtype=np.float)
#     energy = np.linspace(e_min,e_max,nedos)
#     for ieng in range(nedos):
#         for ikpt in range(num_kpoints):
#             for ibnd in range(eigenvals[0].shape[0]):
#                 dos[ieng,0] = energy[ieng]
#                 # print np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
#                 dos[ieng,1] += (coeff[ikpt,ibnd,iwan,0]**2+coeff[ikpt,ibnd,iwan,1]**2)*np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
#     return dos

def get_WOOP(U_matrix, kpoints, R1, R2, iwan1, iwan2, eigenvals, num_kpoints, e_min, e_max, nedos, T=0.1):
    '''

    WOOP = \sum_{j,k} C^*_{iwan1,R1;j} C_{iwan2,R2;j} \delta(E_j-\epsilon)

    returns: dos
    format:  dos[ieng,energy/dos]
    '''
    dos = np.zeros((nedos,2), dtype=float)
    energy = np.linspace(e_min,e_max,nedos)
    for ieng in range(nedos):
        dos[ieng,0] = energy[ieng]
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                dos[ieng,1] += np.real(np.dot(np.conj(coef_gen(U_matrix,ibnd,iwan1,R1,kpoints,ikpt)),coef_gen(U_matrix,ibnd,iwan2,R2,kpoints,ikpt))) \
                               * 1/np.sqrt(np.pi*T)* np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    return dos


def get_WOHP(Hopping, U_matrix, kpoints, R1, R2, iwan1, iwan2, eigenvals, num_kpoints, e_min, e_max, nedos, T=0.2):
    '''
    Hamiltonian weighted projected density of states under Wannier picture.

    WOHP = -H_{iwan1,R1; iwan2,R2} \sum_{j,k} C^*_{iwan1,R1;j} C_{iwan2,R2;j} \delta(E_j-\epsilon)

    returns: dos
    format:  dos[ieng,energy/dos]
    '''
    dos = np.zeros((nedos,2), dtype=float)
    energy = np.linspace(e_min,e_max,nedos)
    for ieng in range(nedos):
        dos[ieng,0] = energy[ieng]
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                dos[ieng,1] += -Hopping*np.real(np.dot(np.conj(coef_gen(U_matrix,ibnd,iwan1,R1,kpoints,ikpt)),coef_gen(U_matrix,ibnd,iwan2,R2,kpoints,ikpt))) \
                               * 1/np.sqrt(np.pi*T)* np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    return dos
