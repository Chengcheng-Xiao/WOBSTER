from __future__ import print_function
import numpy as np

def read_hr(Filename='wannier90_hr.dat'):
    '''
    Read in collinear wannier90_hr.dat and construct a full spinor Hamiltonion by
    filling the upper left (up) or lower right (down) block.
    also return the index of the home cell which we will use to add Hsoc.
    '''
    print('reading '+Filename+ ' ...')
    f=open(Filename,'r')
    data=f.readlines()
    #read hopping matrix
    num_wann = int(data[1])
    nrpts = int(data[2])
    r_hop= np.zeros([num_wann,num_wann], dtype=complex)
    #hop=[]
    #skip n lines of degeneracy of each Wigner-Seitz grid point
    skiplines = int(np.ceil(nrpts / 15.0))
    istart = 3 + skiplines
    deg=[]
    for i in range(3,istart):
        deg.append(np.array([int(j) for j in data[i].split()]))
    deg=np.concatenate(deg,0)

    icount=0
    ii=-1
    Rlatt = []
    hopps = []
    for i in range(istart,len(data)):
        line=data[i].split()
        m = int(line[3]) - 1
        n = int(line[4]) - 1
        r_hop[m,n] = complex(round(float(line[5]),6),round(float(line[6]),6))
        icount+=1
        if(icount % (num_wann*num_wann) == 0):
            ii+=1
            R = np.array([float(x) for x in line[0:3]])
            #hop.append(np.asarray([R,r_hop]))
            #r_hop= np.zeros([num_wann,num_wann], dtype=complex)

            Rlatt.append(R)
            hopps.append(r_hop)
            #hop.append(np.asarray([R,r_hop]))
            r_hop= np.zeros([num_wann,num_wann], dtype=complex)
    Rlatt=np.asarray(Rlatt)
    hopps=np.asarray(hopps)
    deg = np.reshape(deg,[nrpts,1,1])
    hopps=hopps/deg

    hop_fin = np.zeros([nrpts,num_wann,num_wann],dtype=complex)
    for i in range(nrpts):
        shape = np.eye(1)
        hop_fin[i] = np.kron(shape, hopps[i])
    # print('successfully')
    return hop_fin,Rlatt,deg

def read_u_matrix(file_name,dimension_fix,dimension,num_kpoints):

    '''
    a subroutine to get u matrix from seedname_u.mat, optional padding with zeros

    dimension_fix -> total number of bloch bands(can be extended by padding zeros),
    dimension -> number of wannier functions, num_kpoints -> number of kpoints

    returns: matrixs, kpoints
    Format: matrixs[ikpt,iwann,ibnds,real/imaginary]
            kpoints[ikpt,x/y/z]
    '''

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
#                    matrixs[k][num_wann][num_bands][dat]= matrix[dimension*num_bands+num_wann][dat] # Manual tells us "row first", and its a LIE?
                    matrixs[k][num_bands][num_wann][dat]= matrix[dimension*num_wann+num_bands][dat]
    return np.array(matrixs,dtype=np.float),np.array(kpoints,dtype=np.float)

def read_eig(file_name,dimension,num_kpoints):
    '''
    get eigenvalues from wannier90.eig file

    returns: eigenvals
    Format: eigenvals[kpoint_number,band_number]
    '''
    #readin data to data[]
    file = open(file_name,"r")
    content = [x.rstrip("\n") for x in file]
    data = np.array([x.split()[2] for x in content[:]],dtype=np.float)
    eigenvals = data.reshape(num_kpoints,dimension)
    return eigenvals
