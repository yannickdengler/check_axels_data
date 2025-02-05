import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import scattering as result
import math


mpi = 0.38649
mrho = 0.5494

def error_of_array(array):
    tempo = []
    for i in range(len(array[0])):
        tmp = array[:,i]
        tmp.sort()
        tempo.append(tmp)
    num_perc = math.erf(1/np.sqrt(2))
    length = len(array)
    result = np.transpose(tempo)
    return result[length//2], abs(result[length//2]-result[math.floor(length*(1-num_perc)/2)]), abs(result[length//2]-result[math.ceil(length*(1+num_perc)/2)]),result[math.floor(length*(1-num_perc)/2)],result[math.ceil(length*(1+num_perc)/2)]


def en_L(N_L,p2):
    return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
def en_rho_L(N_L,p2):
    return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi

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

def plot_E_prime(file, show=False, save = True, pref = ""):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    res,  res_sample = result.read_from_hdf(file)
    num_gaussian = len(res_sample["E_prime"])

    N_Ls = res["N_Ls"]
    plt.ylabel("E/m$\pi$")
    plt.xlabel("N_L")
    plt.grid()
    xarr = np.linspace(min(N_Ls)-1, max(N_Ls)+1)
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

    E_pipi_prime = error_of_array(res_sample["E_prime"])

    plt.errorbar(N_Ls, E_pipi_prime[0], yerr = [E_pipi_prime[1],E_pipi_prime[2]], capsize=3, ls="", marker = "x")

    plt.legend()
    plt.savefig("en_l_new.pdf")
    plt.show()
    plt.clf()


def plot_E(file, show=False, save = True, pref = ""):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    res,  res_sample = result.read_from_hdf(file)
    num_gaussian = len(res_sample["E_prime"])

    N_Ls = res["N_Ls"]
    plt.ylabel("E/m$\pi$")
    plt.xlabel("N_L")
    plt.grid()
    xarr = np.linspace(min(N_Ls)-1, max(N_Ls)+1)
    pipi1arr = [mpi*en_L(x,1) for x in xarr]
    pipi2arr = [mpi*en_L(x,2) for x in xarr]
    rho1arr = [mpi*en_rho_L(x,1) for x in xarr]
    rho2arr = [mpi*en_rho_L(x,2) for x in xarr]
    plt.axhline(mpi*2, color = "green", label = "pipi")
    plt.axhline(mrho, color = "orange", label = "rho")
    plt.plot(xarr, pipi1arr, ls = "dashdot", color = "green")
    plt.plot(xarr, pipi2arr, ls = "dashdot", color = "green")
    plt.plot(xarr, rho1arr, ls = "dashdot", color = "orange")
    plt.plot(xarr, rho2arr, ls = "dashdot", color = "orange")

    E_pipi_prime = [mpi*x for x in error_of_array(res_sample["E_prime"])]

    plt.errorbar(N_Ls, E_pipi_prime[0], yerr = [E_pipi_prime[1],E_pipi_prime[2]], capsize=3, ls="", marker = "x")

    plt.legend()
    plt.savefig("aen_l_new.pdf")
    plt.show()
    plt.clf()
    
def sqrt_s(E,P):    
    return np.sqrt(E**2-P**2)
def sqrt_s_pi_L(N_L,p2):
    return sqrt_s(1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2),np.sqrt(p2)*2*np.pi/(N_L*mpi))

def plot_sqrt_s(file, show=False, save = True, pref = ""):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    res,  res_sample = result.read_from_hdf(file)
    num_gaussian = len(res_sample["E_prime"])

    N_Ls = res["N_Ls"]    
    plt.ylabel("$\sqrt{s}$/$m_\pi$")
    plt.xlabel("$N_L$")
    plt.grid()
    xarr = np.linspace(min(N_Ls)-1, max(N_Ls)+1)
    pipi1arr = [sqrt_s_pi_L(x,1) for x in xarr]
    pipi2arr = [sqrt_s_pi_L(x,2) for x in xarr]
    plt.axhline(2, color = "green", label = "pipi")
    plt.axhline(mrho/mpi, color = "orange", label = "rho")
    plt.plot(xarr, pipi1arr, ls = "dashdot", color = "green")
    plt.plot(xarr, pipi2arr, ls = "dashdot", color = "green")

    E_cm_prime = error_of_array(res_sample["E_cm_prime"])

    plt.errorbar(N_Ls, E_cm_prime[0], yerr = [E_cm_prime[1],E_cm_prime[2]], capsize=3, ls="", marker = "x")

    plt.legend()
    plt.savefig("sqrt_s_new.pdf")
    plt.show()
    plt.clf()

colors = ["magenta", "red", "green", "black", "orange", "cyan", "olive", "gray", "purple"]

def delete_steps(arr):
    for i in range(len(arr)-1):
        print(arr[i], arr[i+1])
        if arr[i+1] < arr[i]: 
            arr[i] = np.nan
    return arr

def p3_cot_PS(file, show=False, save = True, pref = "", x_ax = "sqrt_s"):
    plt.rcParams['figure.figsize'] = [10, 6]
    fontsize = 14
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    fig, ax = plt.subplots()
    res,  res_sample = result.read_from_hdf(file)
    num_gaussian = len(res_sample["E_prime"])

    N_Ls = res["N_Ls"]    
    plt.ylabel("$p^3\, \cot(\delta)/m_\pi^3$")
    plt.grid()

    if x_ax == "sqrt_s":
        plt.xlabel("$\sqrt{s}/m_\pi$")
        E_cm_prime = res["E_cm_prime"]
        E_cm_prime_sam = np.transpose(res_sample["E_cm_prime"])
        ax.set_xlim([2,3.2])
    elif x_ax == "aE":
        plt.xlabel("aE")
        E_cm_prime = res["E_pure"]
        E_cm_prime_sam = np.transpose(res_sample["E_pure"])
        ax.set_xlim([0.8, 1.3])


    P3_cot_PS_prime = np.real(res["P3_cot_PS_prime"])
    P3_cot_PS_prime_sam = np.transpose(np.real(res_sample["P3_cot_PS_prime"]))

    length = len(E_cm_prime_sam[0])

    lvs = [1,2,1,2,1,2,1,1,2]
    Ps = [1,1,1,1,1,1,2,2,2]

    num_perc = math.erf(1/np.sqrt(2))
    # for i in range(len(E_cm_prime_sam)):
    for i in [1,3,5,8]:
        print(i, E_cm_prime[i])
        ax.scatter(E_cm_prime[i],P3_cot_PS_prime[i], color = colors[i], label = "Lv: %i, |P|=%i, NL=%i"%(lvs[i],Ps[i],N_Ls[i]))
        sorted_indices = np.argsort(E_cm_prime_sam[i])  
        ax.plot(E_cm_prime_sam[i][sorted_indices][math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)],delete_steps(P3_cot_PS_prime_sam[i][sorted_indices])[math.floor(length*(1-num_perc)/2):math.ceil(length*(1+num_perc)/2)], color = colors[i])


    ax.set_ylim([-5,5])
    ax.legend(loc="best")
    fig.savefig("p3_cotPS_%s.pdf"%(x_ax))
    plt.show()
    fig.clf()


if __name__ == "__main__":
    # plot_E("PS_69_092")    
    # plot_sqrt_s("PS_69_092")   
    p3_cot_PS("PS_69_092")    
    p3_cot_PS("PS_69_092", x_ax="aE")    

