import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm


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
    return 2*np.sqrt(mpi**2+(q_of_E(E, dtot, N_L, mpi)*(2*np.pi)/N_L)**2)/mpi

def phase_shift(Lambda, mpi, E, dtot, N_L):
    return (np.arctan(1/cot_delta_Lambda(norm2(dtot),Lambda)(mpi,mpi,q_of_E(E,dtot,N_L,mpi)**2,N_L).real)*360/(2*np.pi))%180

def P_vecf(i):
    if i == 1:
        return [0,0,1]
    else:
        return [1,1,0]

data = np.transpose(np.genfromtxt("data/axel_data.dat"))



print(data.shape)

mpi = 0.38649
mrho = 0.5494

NL_arr = [int(x) for x in data[0]]
P_arr = data[4]
en_lv_arr = data[5]
en_arr = data[6]/mpi
en_p_arr = data[7]/mpi
en_n_arr = data[8]/mpi

print(NL_arr)

en_sample_arr = []

sarr = []
PcotPSarr = []

N_resample = 10

for i in range(len(NL_arr)):
    if en_arr[i] > 2*mpi:
        en_sample_arr.append([])
        sarr.append([])
        PcotPSarr.append([])
        for E in np.linspace(en_arr[i]-en_n_arr[i],en_arr[i],N_resample):
            en_sample_arr[i].append(E)
        for E in np.linspace(en_arr[i],en_arr[i]+en_p_arr[i],N_resample):
            # en_sample_arr[i].append(max(2*mpi*1.01,E))
            en_sample_arr[i].append(E)
    if i < 5.5:
        for j in range(len(en_sample_arr[i])):
            sarr[i].append(s_of_E(en_sample_arr[i][j],P_vecf(P_arr[i]),NL_arr[i],mpi))
            if sarr[i][j] > 2:
                PcotPSarr[i].append(cot_delta_001(mpi,mpi,q_of_E(en_sample_arr[i][j],P_vecf(P_arr[i]),NL_arr[i],mpi)**2,NL_arr[i]))
            else:
                PcotPSarr[i].append(0)
    else: 
        for j in range(len(en_sample_arr[i])):
            sarr[i].append(s_of_E(en_sample_arr[i][j],P_vecf(P_arr[i]),NL_arr[i],mpi))
            if sarr[i][j] > 2:
                print(en_sample_arr[i][j],P_vecf(P_arr[i]),NL_arr[i],mpi)
                PcotPSarr[i].append(cot_delta_110(mpi,mpi,q_of_E(en_sample_arr[i][j],P_vecf(P_arr[i]),NL_arr[i],mpi)**2,NL_arr[i]))
            else:
                PcotPSarr[i].append(0)

plt.xlabel("Ecm/mpi")
plt.ylabel("cotPS")

sarrplot = np.real(sarr)
PcotPSarrplot = np.real(PcotPSarr)


for i in range(len(NL_arr)):
    plt.plot(sarrplot[i], PcotPSarrplot[i], label = "L%i,P%i,lv%i"%(NL_arr[i],P_arr[i],en_lv_arr[i]))
plt.legend()
plt.savefig("PS_BS.pdf")
plt.show()
        







# plt.grid()
# xarr = np.linspace(min(NL_arr)-1, max(NL_arr)+1)
# pipi1arr = [en_L(x,1) for x in xarr]
# pipi2arr = [en_L(x,2) for x in xarr]
# rho1arr = [en_rho_L(x,1) for x in xarr]
# rho2arr = [en_rho_L(x,2) for x in xarr]
# plt.axhline(2, color = "green", label = "pipi")
# plt.axhline(mrho/mpi, color = "orange", label = "rho")
# plt.plot(xarr, pipi1arr, ls = "dashdot", color = "green")
# plt.plot(xarr, pipi2arr, ls = "dashdot", color = "green")
# plt.plot(xarr, rho1arr, ls = "dashdot", color = "orange")
# plt.plot(xarr, rho2arr, ls = "dashdot", color = "orange")

# def marker(pind):
#     if pind == 1:
#         return "o"
#     else:
#         return "x"
# def colorf(en_lv_ind):
#     if en_lv_ind == 1:
#         return "blue"
#     else:
#         return "red"

# for i in range(len(en_arr)):
#     plt.errorbar(NL_arr[i], en_arr[i], yerr = [[en_n_arr[i],],[en_p_arr[i],]],marker = marker(P_arr[i]), color = colorf(en_lv_arr[i]), capsize=3)

# plt.legend()
# plt.show()