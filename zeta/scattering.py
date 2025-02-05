from pylink_wlm import wlm_func_c as wlm
import numpy as np

import matplotlib.pyplot as plt

I = complex(0,1)

# Definitions from arXiv:1105.5636v3 or arXiv:1011.5288v2 or arXiv:2012.09761v2


def cot_delta_001(m1,m2,q2,N_L):
    return wlm(0,0,0,0,1,m1,m2,q2,N_L) + 2*wlm(2,0,0,0,1,m1,m2,q2,N_L)
def cot_delta_110(m1,m2,q2,N_L):
    return wlm(0,0,1,1,0,m1,m2,q2,N_L) - wlm(2,0,1,1,0,m1,m2,q2,N_L) + np.sqrt(3/2) * (wlm(2,2,1,1,0,m1,m2,q2,N_L) - wlm(2,-2,1,1,0,m1,m2,q2,N_L))

def norm(vec):
    return np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
def norm2(vec):
    return vec[0]**2+vec[1]**2+vec[2]**2

def q_of_E(E, dtot, N_L, mpi):
    return np.sqrt((E**2-(norm(dtot)*2*np.pi/N_L)**2)/4-mpi**2)*N_L/(2*np.pi)

def s_of_E(E, dtot, N_L, mpi):
    return 2*np.sqrt(mpi**2+(q_of_E(E, dtot, N_L, mpi)*(2*np.pi)/N_L)**2)

def phase_shift(Lambda, mpi, E, dtot, N_L):
    return (np.arctan(1/cot_delta_Lambda(norm2(dtot),Lambda)(mpi,mpi,q_of_E(E,dtot,N_L,mpi)**2,N_L).real)*360/(2*np.pi))%180

if __name__ == "__main__":

    print("b6.9_m-0.92")

    mpi = 0.38649
    
    data = np.genfromtxt("data/axel_data.dat")

    data_b69m092 = data[:5]
    energy_levels_b69_m092 = data_b69m092[::,6:]


    NL_arr = [16,16,14,24,24]
    P_arr = [1,2,1,1,2]
    Pvec_arr = [[0,0,1],[1,1,0],[0,0,1],[0,0,1],[1,1,0]]

    cot_PS_arr = []
    s_arr = []

    for i in range(len(energy_levels_b69_m092)):
        print(i)
        if energy_levels_b69_m092[i,3] != 0:
            if P_arr[i] == 1:
                # s_arr.append(s_of_E(energy_levels_b69_m092[i,0],Pvec_arr[i],NL_arr[i],mpi))
                # cot_PS_arr.append(cot_delta_001(mpi,mpi,q_of_E(energy_levels_b69_m092[i,0],Pvec_arr[i],NL_arr[i],mpi)**2,NL_arr[i]))
                s_arr.append(s_of_E(energy_levels_b69_m092[i,3],Pvec_arr[i],NL_arr[i],mpi))
                cot_PS_arr.append(cot_delta_001(mpi,mpi,q_of_E(energy_levels_b69_m092[i,3],Pvec_arr[i],NL_arr[i],mpi)**2,NL_arr[i]))
            elif P_arr[i] == 2:
                
                # s_arr.append(s_of_E(energy_levels_b69_m092[i,0],Pvec_arr[i],NL_arr[i],mpi))
                # cot_PS_arr.append(cot_delta_110(mpi,mpi,q_of_E(energy_levels_b69_m092[i,0],Pvec_arr[i],NL_arr[i],mpi)**2,NL_arr[i]))
                s_arr.append(s_of_E(energy_levels_b69_m092[i,3],Pvec_arr[i],NL_arr[i],mpi))
                cot_PS_arr.append(cot_delta_110(mpi,mpi,q_of_E(energy_levels_b69_m092[i,3],Pvec_arr[i],NL_arr[i],mpi)**2,NL_arr[i]))
            else:
                print("wrong momenutm")
                exit()
    
    print(s_arr)
    print(cot_PS_arr)

    plt.scatter(s_arr, cot_PS_arr)
    plt.show()