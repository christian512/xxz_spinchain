import numpy as np
from scipy.sparse import csr_matrix

def runge_kutta_vec(mat, vector,h,t_end,t_start=0):

    steps = int((t_end-t_start)/h)
    ret_list = [None]*(steps+1)
    t_list = np.empty(steps+1)
    ret_list[0] = vector
    t_list[0] = t_start
    t = t_start
    for i in range(steps):
        #print('runge_kutta: ' + str(i) + ' / ' + str(steps))
        k1 = np.dot(mat,vector)
        k2 = np.dot(mat,vector+h/2*k1)
        k3 = np.dot(mat,vector+h/2*k2)
        k4 = np.dot(mat,vector+h*k3)
        vector = vector + h/6*(k1+2*k2+2*k3+k4)
        t = t + h

        #check norm of vector
        #print(np.linalg.norm(vector))
        if np.abs(np.linalg.norm(vector)-1) > 1e-2:
            print('Norm is not conserved, restarting Runge-Kutta with h= ' + str(h/2))
            return runge_kutta_vec(mat,ret_list[0],h/2,t_end,t_start)

        ret_list[i+1] = vector
        t_list[i+1] = t
    return (t_list,ret_list)
#Mat needs to be csr_matrix not np.ndarray
def runge_kutta_vec_onerun(mat,vector,h,t_start):
    k1 = mat.dot(vector)
    k2 = mat.dot(vector + h / 2 * k1)
    k3 = mat.dot(vector + h / 2 * k2)
    k4 = mat.dot(vector + h * k3)
    vector = vector + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    t = t_start + h

    if np.abs(np.linalg.norm(vector) - 1) > 1e-2:
        print('Norm is not conserved!')
        return None
    return (t,vector)
