import numpy as np
from np_helper import state_to_num,indexOfState
from scipy.sparse import lil_matrix


states_num_glob = []

"""Constructing the hamiltonian for the reduced size
    Input:  states_np   : np.array  : States in numpy array form
            J           : float     : Coupling constant for spin flips
            r           : float     : Coupling constant for s_z term
            num_part    : num_part  : Number of particles
    Output: ham         : csr-arr   : Hamiltonian of the reduced size in a sparse array format   
"""
def constructHamilton(states_num,J,r,epsilons,num_part):
    #save the states to global variable
    global states_num_glob
    states_num_glob = states_num



    #define the hamiltonian
    ham = lil_matrix((states_num_glob.shape[0],states_num_glob.shape[0]))

    #loop through all states
    for i in range(states_num_glob.shape[0]):
        state_num = states_num_glob[i]

        #place defects
        for j in range(num_part):
            ham[i,i] += epsilons[j] * sz(num_part-j,state_num)

        for j in range(num_part-1):
            ham[i,i] -= r*sz_sz(j,j+1,state_num)
            res = sp_sm(j,j+1,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J/2
            res = sm_sp(j,j+1,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J/2
    return ham.tocsr()

def constructHamiltonPeriodic(states_num,J,r,epsilons,num_part):
    #save the states to global variable
    global states_num_glob
    states_num_glob = states_num



    #define the hamiltonian
    ham = lil_matrix((states_num_glob.shape[0],states_num_glob.shape[0]))

    #loop through all states
    for i in range(states_num_glob.shape[0]):
        state_num = states_num_glob[i]

        #place defects
        for j in range(num_part):
            ham[i,i] += J*epsilons[j] * sz(num_part-j,state_num)

        for j in range(num_part):
            l = j + 1
            if l == num_part:
                l = 0
            ham[i,i] -= J*r*sz_sz(j,l,state_num)
            res = sp_sm(j,l,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J
            res = sm_sp(j,l,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J
    return ham.tocsr()

def constructHamiltonNNN(states_num,J,r,J2,r2,epsilons,num_part):
    #First construct NN hamiltonian
    # save the states to global variable
    global states_num_glob
    states_num_glob = states_num

    # define the hamiltonian
    ham = lil_matrix((states_num_glob.shape[0], states_num_glob.shape[0]))

    # loop through all states
    for i in range(states_num_glob.shape[0]):
        state_num = states_num_glob[i]

        # place defects
        for j in range(num_part):
            ham[i, i] +=  epsilons[j] * sz(num_part - j, state_num)
        # NN interaction
        for j in range(num_part - 1):
            ham[i, i] -= r * sz_sz(j, j + 1, state_num)
            res = sp_sm(j, j + 1, state_num)
            if not type(res) == type(None):
                k = indexOfState(res, states_num_glob)[0][0]
                ham[k, i] -= J/2
            res = sm_sp(j, j + 1, state_num)
            if not type(res) == type(None):
                k = indexOfState(res, states_num_glob)[0][0]
                ham[k, i] -= J/2

        #NNN interaction
        for j in range(num_part-2):
            ham[i,i] -= r2*sz_sz(j,j+2,state_num)
            res = sp_sm(j,j+2,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J2/2
            res = sm_sp(j,j+2,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J2/2
    return ham.tocsr()

def constructHamiltonPeriodicNNN(states_num,J,r,J2,r2,epsilons,num_part):
    #First construct NN hamiltonian
    # save the states to global variable
    global states_num_glob
    states_num_glob = states_num

    # define the hamiltonian
    ham = lil_matrix((states_num_glob.shape[0], states_num_glob.shape[0]))

    # loop through all states
    for i in range(states_num_glob.shape[0]):
        state_num = states_num_glob[i]

        # place defects
        for j in range(num_part):
            ham[i, i] += J * epsilons[j] * sz(num_part - j, state_num)
        # NN interaction
        for j in range(num_part):
            l = j + 1
            if l == num_part:
                l = 0
            ham[i, i] -= J * r * sz_sz(j, l, state_num)
            res = sp_sm(j, l, state_num)
            if not type(res) == type(None):
                k = indexOfState(res, states_num_glob)[0][0]
                ham[k, i] -= J
            res = sm_sp(j, l, state_num)
            if not type(res) == type(None):
                k = indexOfState(res, states_num_glob)[0][0]
                ham[k, i] -= J

        #NNN interaction
        for j in range(num_part):
            l = j + 2
            if j == num_part-2:
                l = 0
            if j == num_part-1:
                l = 1
            ham[i,i] -= J2*r2*sz_sz(j,l,state_num)
            res = sp_sm(j,l,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J2
            res = sm_sp(j,l,state_num)
            if not type(res) == type(None):
                k = indexOfState(res,states_num_glob)[0][0]
                ham[k,i] -= J2
    return ham.tocsr()

"""Constructing the magnetization matrix for given states and given particle list"""
def constructMagnetization(states_num,particles):
    #define matrix
    mag = lil_matrix((states_num.shape[0],states_num.shape[0]))
    #set elements
    for i in range(states_num.shape[0]):
        for j in particles:
            mag[i,i] += sz(j,states_num[i])
    return mag.tocsr()

    return 0

""" Function that calculates the spin operator s_z and s_z on particles i and j 
    to a given state in decimal number representation
    WATCH OUT: i and j are counted from behind"""
def sz_sz(i,j,state_num):
    #if you give a list of states to the function
    if type(state_num) == np.ndarray:
        res = np.empty(state_num.shape[0])
        for k in range(res.shape[0]):
            res[k] = sz_sz(i,j,state_num[k])
        return res
    #if both are spin up
    if np.bitwise_and(2**i,state_num) == 2**i and np.bitwise_and(2**j,state_num) == 2**j: return 1/4
    elif np.bitwise_and(2**i,state_num) == 0 and np.bitwise_and(2**j,state_num) == 0: return 1/4
    else: return -1/4

""" Function that calculates the spin operator S_+ and S_- on particles i and j 
    to a given state in decimal number representation
    WATCH OUT: i and j are counted from behind"""
def sp_sm(i,j,state_num):
    #spin-i down and spin-j up
    if np.bitwise_and(2**i,state_num) == 0 and np.bitwise_and(2**j,state_num) == 2**j:
        return np.bitwise_xor(2**i+2**j,state_num)
    return None


""" Function that calculates the spin operator S_- and S_+ on particles i and j 
    to a given state in decimal number representation
    WATCH OUT: i and j are counted from behind"""
def sm_sp(i,j,state_num):
    #spin-i up and spin-j down
    if np.bitwise_and(2**i,state_num) == 2**i and np.bitwise_and(2**j,state_num) == 0:
        return np.bitwise_xor(2**i+2**j,state_num)
    return None


""" Function that calculates the spin operator S_z on particle i to a given state in decimal number representation
    WATCH OUT: i is counted from behind"""
def sz(i,state_num):
    if np.bitwise_and(2**i,state_num) == 2**i: return 1/2
    else: return -1/2


"""Calculating the partial trace on a given coefficient array
    coeffs is a 2d numpy array which contains all coeffs for different times
"""
def partial_trace(dens_mat, j_list,num_parts,states):
    states_num_glob = states
    if type(j_list) == int:
        j_list = np.array([j_list])
    #copy the list of all states
    new_states_num = np.copy(states_num_glob)
    #create a list that contains all deleted parts
    del_states = np.empty(states_num_glob.shape[0])
    #run through all states
    for i in range(states_num_glob.shape[0]):
        #print('i=' + str(i))
        #deleted part
        deleted_part = ''
        s = np.binary_repr(new_states_num[i])
        #add starting zeros
        if len(s) < num_parts:
            s = '0'*(num_parts-len(s))+s
        #run through all particles which should be traced out
        for k in range(j_list.shape[0]):
            #delete the part and store the deleted part and new state in corresponding arrays
            deleted_part += s[j_list[k]-k]
            #print('after s=' + str(s))
            s = s[:j_list[k]-k] + s[(j_list[k]-k+1):]
            #print('after s=' + str(s))
        new_states_num[i] = int(s,2)
        del_states[i] = int(deleted_part,2)

    #create a copy of the new_states to store where the index where the old states gone to
    new_indices = np.copy(new_states_num)

    #delete all duplicates from the new states
    new_states_num = np.unique(new_states_num)

    #write indices of new states to the indice array
    for i in range(new_indices.shape[0]):
        new_indices[i] = np.where(new_states_num==new_indices[i])[0][0]

    #create new reduced density matrix
    red_dens = np.zeros([new_states_num.shape[0],new_states_num.shape[0]],dtype=complex)

    for i in range(states_num_glob.shape[0]):
        for j in range(states_num_glob.shape[0]):
            if del_states[i] == del_states[j]:
                #print('found same deleted: ' + str(i) + '   ' + str(j))
                k = new_indices[i]
                l = new_indices[j]
                red_dens[k,l] += dens_mat[i,j]
                #print(new_indices[i])
    return red_dens

def density_matrix(coeffs):
    return np.outer(coeffs,np.conjugate(coeffs))

def entropy_vn(rho):
    eigvals = np.linalg.eigvals(rho)
    eigvals = eigvals[eigvals > 0]
    sum_list = - eigvals*np.log(eigvals)
    return np.sum(sum_list)

"""Calculates the trace distance between two densitiy matrices"""
def trace_dist(rho1,rho2):
    rho = rho1-rho2
    rho = np.sqrt(np.dot(np.conjugate(np.transpose(rho)),rho))
    eigvals = np.linalg.eigvals(rho)
    return 1/2*np.sum(eigvals)
