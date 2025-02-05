import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from scipy.optimize import bisect

# xarr = np.linspace13,25]

def cot_delta_001(m1,m2,q2,N_L):
    return wlm(0,0,0,0,1,m1,m2,q2,N_L) + 2*wlm(2,0,0,0,1,m1,m2,q2,N_L)
def cot_delta_110(m1,m2,q2,N_L):
    return wlm(0,0,1,1,0,m1,m2,q2,N_L) - wlm(2,0,1,1,0,m1,m2,q2,N_L) + np.sqrt(3/2) * (wlm(2,2,1,1,0,m1,m2,q2,N_L) - wlm(2,-2,1,1,0,m1,m2,q2,N_L))

data = np.transpose(np.genfromtxt("data/axel_data.dat"))

print(data.shape)

mpi = 0.38649
mrho = 0.5494

NL_arr = data[0]
P_arr = data[4]
en_lv_arr = data[5]
en_arr = data[6]/mpi
en_p_arr = data[7]/mpi
en_n_arr = data[8]/mpi

def ecm(e,P):
    return np.sqrt(e**2+P*+2)

xarr = np.linspace(0.5,1.5)
yarr = [cot_delta_001(mpi,mpi,x,24) for x in xarr]

plt.plot(xarr,yarr)
plt.show()

def pcot001_zero(q2,NL):
    return cot_delta_001(mpi,mpi,q2,N_L)

xarr = np.linspace(min(NL_arr)-1, max(NL_arr)+1)

yarr = [bisect(pcot001_zero,0.2,1.8,args=(x,)) for x in xarr]

plt.plt(xarr,yarr)


















# m1=0.5
# m2=m1

# N_L = 16

# xarr = np.linspace(0.5,1.5)
# yarr0 = [np.real(wlm(0,0,0,0,0,m1,m2,x,N_L)) for x in xarr]
# yarr1 = [np.real(wlm(0,0,0,0,1,m1,m2,x,N_L)) for x in xarr]
# yarr2 = [np.real(wlm(2,0,0,0,1,m1,m2,x,N_L)) for x in xarr]
# yarr3 = [np.real(wlm(2,0,1,1,0,m1,m2,x,N_L)) for x in xarr]
# yarr4 = [np.imag(wlm(2,2,1,1,0,m1,m2,x,N_L)) for x in xarr]
# yarr5 = [np.imag(wlm(2,-2,1,1,0,m1,m2,x,N_L)) for x in xarr]

# plt.plot(xarr, yarr0, label="00 000")
# print(yarr0)
# plt.plot(xarr, yarr1, label="00 001")
# print(yarr1)
# plt.plot(xarr, yarr2, label="20 001")
# print(yarr2)
# plt.plot(xarr, yarr3, label="20 110")
# print(yarr3)
# plt.plot(xarr, yarr4, label="22 110")
# print(yarr4)
# plt.plot(xarr, yarr5, label="2m2 110")
# print(yarr5)
# plt.legend()
# plt.show()
# exit()





def en_L(N_L,p2):
    return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
def en_rho_L(N_L,p2):
    return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi

plt.grid()
xarr = np.linspace(min(NL_arr)-1, max(NL_arr)+1)
pipi1arr = [en_L(x,1) for x in xarr]
pipi2arr = [en_L(x,2) for x in xarr]
rho1arr = [en_rho_L(x,1) for x in xarr]
rho2arr = [en_rho_L(x,2) for x in xarr]
plt.axhline(2, color = "green", label = "pipi")
plt.axhline(mrho/mpi, color = "orange", label = "rho")
plt.plot(xarr, pipi1arr, ls = "dashdot", color = "green")
plt.plot(xarr, pipi2arr, ls = "dashdot", color = "green")
plt.plot(xarr, rho1arr, ls = "dashdot", color = "orange")
plt.plot(xarr, rho2arr, ls = "dashdot", color = "orange")

def wlm_zero(q2, N_L):
    return wlm(0,0,0,0,1,mpi,mpi,q2,N_L) + 2*wlm(2,0,0,0,1,m1,m2,q2,N_L)

zeta_zero = [bisect(wlm_zero,1e-5,1-1e-5, args=L) for L in xarr]

def marker(pind):
    if pind == 1:
        return "o"
    else:
        return "x"
def colorf(en_lv_ind):
    if en_lv_ind == 1:
        return "blue"
    else:
        return "red"

for i in range(len(en_arr)):
    plt.errorbar(NL_arr[i], en_arr[i], yerr = [[en_n_arr[i],],[en_p_arr[i],]],marker = marker(P_arr[i]), color = colorf(en_lv_arr[i]), capsize=3)

plt.title("red exited, 001-o, 110-x")

plt.legend()
plt.savefig("en_l.pdf")
plt.show()