import numpy as np
import matplotlib.pyplot as plt

def get_u_matrix(file_name,dimension_fix,dimension,num_kpoints):

    ''' a subroutine to get u matrix from seedname_u.mat, optional padding provided '''
    ''' dimension_fix -> total number of bloch bands(can be extended by padding zeros), dimension -> number of wannier functions, num_kpoints -> number of kpoints '''

    #readin data to data[]
    file = open(file_name,"r")
    content = [x.rstrip("\n") for x in file]
    data = [x.split()[:4] for x in content[3:]]
    #initialise array
    matrix = np.zeros((dimension_fix**2*num_kpoints,2)) #initialize an empty array with zeros.
    matrixs  = np.zeros((num_kpoints,dimension_fix,dimension,2))
    kpoints = np.zeros((num_kpoints,3))
    #labeling each dataset
    #######################################
    #the format of seedname_u.mat is
    # U_band1_wann1[rel] U_band1_wann1[img]
    # U_band2_wann1[rel] U_band2_wann1[img]
    # U_band3_wann1[rel] U_band3_wann1[img]
    # ...
    # U_bandn_wannm[rel] U_bandn_wannm[img]
    #######################################
    #reading Kpoints
    for k in range(num_kpoints):
        kpoints[k]=data[k*(dimension*dimension+2)]
    #differentiating img and rel part of data
    for k in range(num_kpoints):
        for i in range(k*(dimension**2+2)+1,k*(dimension**2+2)+(dimension**2+1)):
            for dat in range(2):
                matrix[i-(k*(dimension**2+2))-1][dat]=data[i][dat]
    #labeling bloch band and wannier band number
        for num_wann in range(dimension):
            for num_bands in range(dimension):
                for dat in range(2):
#                    matrixs[k][num_wann][num_bands][dat]= matrix[dimension*num_bands+num_wann][dat] # Manual tells us "row first", and its a LIE
                    matrixs[k][num_bands][num_wann][dat]= matrix[dimension*num_wann+num_bands][dat]
    return matrixs,kpoints

def get_eig(file_name,dimension,num_kpoints):
    '''
    get eigenvalues from wannier90.eig file
    '''
    #readin data to data[]
    file = open(file_name,"r")
    content = [x.rstrip("\n") for x in file]
    data = np.array([x.split()[2] for x in content[:]],dtype=np.float)
    # the format is eigenvals[kpoint_number,band_number]
    eigenvals = data.reshape(num_kpoints,dimension)
    return eigenvals

def get_dos(eigenvals, num_kpoints, e_min, e_max, nedos, T=0.1):
    '''
    calculate density of states by:

    DOS(\epsilon) = \sum_j \delta(E_j-\epsilon)

    \delta function is expressed as:
    \delta = exp(-(E_j-\epsilon)**2/T)

    where T is the semaring factor.
    '''
    dos = np.zeros((nedos,2), dtype=np.float)
    energy = np.linspace(e_min,e_max,nedos)
    for ieng in range(nedos):
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                dos[ieng,0] = energy[ieng]
                # print np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
                dos[ieng,1] += np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    return dos

def coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt):
    '''
    U_matrix: full unitary matrix
    ibnd: band number
    iwan: wannier function number
    R: lattice vectore, should be array(3)
    kpoints: array contains all kpoints info
    ikpt: number kpoint
    '''
    coef = np.conj(complex(U_matrix[ikpt,ibnd,iwan,0],U_matrix[ikpt,ibnd,iwan,1]))
    kdR = 2*np.pi*np.dot(kpoints[ikpt],R)
    #coef *= ((2*np.pi)**3/V)*np.exp(complex(0,kdR))
    coef *= np.exp(complex(0,kdR))
    return coef

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

def get_WOOP(U_matrix, kpoints, R, iwan, eigenvals, num_kpoints, e_min, e_max, nedos, T=0.1):
    '''
    WOOP under Wannier picture = projected DOS

    WO)P = \sum_{j,k} C^*_{iwan,R;j} C_{iwan,R;j} \delta(E_j-\epsilon)

    '''
    dos = np.zeros((nedos,2), dtype=np.float)
    energy = np.linspace(e_min,e_max,nedos)
    # print 'fuck'
    for ieng in range(nedos):
        dos[ieng,0] = energy[ieng]
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                # if ieng == 0 and ikpt ==0:
                #     print ieng,ikpt,ibnd
                #     print np.real(np.dot(np.conj(coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt)),coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt)))
                dos[ieng,1] += np.real(np.dot(np.conj(coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt)),coef_gen(U_matrix,ibnd,iwan,R,kpoints,ikpt))) \
                               * np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
                # dos[ieng,1] += (coeff[ikpt,ibnd,iwan,0]**2+coeff[ikpt,ibnd,iwan,1]**2)*np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    # print 'fuck'
    return dos


def get_WOHP(Hopping, U_matrix, kpoints, R1, R2, iwan1, iwan2, eigenvals, num_kpoints, e_min, e_max, nedos, T=0.2):
    '''
    Hamiltonian weighted projected density of states under Wannier picture.

    WOHP = -H_{iwan1,R1; iwan2,R2} \sum_{j,k} C^*_{iwan1,R1;j} C_{iwan2,R2;j} \delta(E_j-\epsilon)
    '''
    dos = np.zeros((nedos,2), dtype=np.float)
    energy = np.linspace(e_min,e_max,nedos)
    # print 'fuck'
    for ieng in range(nedos):
        dos[ieng,0] = energy[ieng]
        for ikpt in range(num_kpoints):
            for ibnd in range(eigenvals[0].shape[0]):
                # print ieng,ikpt,ibnd
                # print np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
                dos[ieng,1] += -Hopping*np.real(np.dot(np.conj(coef_gen(U_matrix,ibnd,iwan1,R1,kpoints,ikpt)),coef_gen(U_matrix,ibnd,iwan2,R2,kpoints,ikpt))) \
                               *np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
                # dos[ieng,1] += (coeff[ikpt,ibnd,iwan,0]**2+coeff[ikpt,ibnd,iwan,1]**2)*np.exp(-(energy[ieng]-eigenvals[ikpt,ibnd])**2/T)/num_kpoints
    # print 'fuck'
    return dos
