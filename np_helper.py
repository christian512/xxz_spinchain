import numpy as np


def state_to_num(state):
    """Converts a state in a binary array form [0,1,0,..,1] to a integer"""
    if len(state.shape) == 2:
        c = 2 ** np.arange(state.shape[1])[::-1]
        nums = np.empty(state.shape[0],dtype=int)
        for j in range(state.shape[0]):
            nums[j] = int(state[j].dot(c))
        return nums
    c = 2**np.arange(state.shape[0])[::-1]
    num = state.dot(c)
    return int(num)


def indexOfState(state_num,states):
    """Finds the index of a specific state in the states array"""
    if states == []:
        print('no states provided ')
        return None
    index = np.where(states == state_num)
    return index



def store_states(states,num_part,filename):
    """Storing the states to a file. Also needs number of particles num_part."""
    f = open(filename,'w')
    for i in range(states.shape[0]):
        f.write(str(i) + '\t :' + bin_rep(states[i],num_part))
        f.write('\n')
    f.close()
    return True


def bin_rep(state, num_part):
    """"Corrects binary representation of a state given the number of particles"""
    s = np.binary_repr(state)
    if len(s) < num_part:
        s = '0'*(num_part-len(s))+ s
    return s


def read_states(file):
    """Reads states from file"""
    f = open(file, 'r')
    content = f.readlines()
    for i in range(len(content)):
        content[i] = content[i].split(':')[1]
        content[i] = content[i][:-1]
        content[i] = list(content[i])
        content[i] = list(map(int,content[i]))
    states = state_to_num(np.array(content))
    return states
