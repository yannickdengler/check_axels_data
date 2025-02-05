import numpy as np
import matplotlib.pyplot as plt

data = np.transpose(np.genfromtxt("data/axel_data.dat"))

print(data.shape)

mpi = 0.38649
mrho = 0.5494

NL_arr = data[0]
P_arr = data[4]
Pvec_arr = np.transpose(data[1:4])
P_m_arr = [np.sqrt(np.dot(Pvec_arr[i],Pvec_arr[i]))*2*np.pi/(mpi*NL_arr[i]) for i in range(len(Pvec_arr))]
en_lv_arr = data[5]
en_arr = data[6]/mpi
en_m_arr = data[7]/mpi
en_p_arr = data[8]/mpi

# print(P_m_arr)
# exit()


# def sqrt_s(E,Pvec):    
#     return np.sqrt(E**2-Pvec[0]*Pvec[0]-Pvec[1]*Pvec[1]-Pvec[2]*Pvec[2])
def sqrt_s(E,P):    
    return np.sqrt(E**2-P**2)

sqrt_s_arr = [sqrt_s(en_arr[i],P_m_arr[i]) for i in range(len(en_arr))]
sqrt_s_m_arr = [sqrt_s(en_arr[i]-en_m_arr[i],P_m_arr[i]) for i in range(len(en_arr))]
sqrt_s_p_arr = [sqrt_s(en_arr[i]+en_p_arr[i],P_m_arr[i]) for i in range(len(en_arr))]


def plot_sqrt_s():
    def sqrt_s_pi_L(N_L,p2):
        return sqrt_s(1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2),np.sqrt(p2)*2*np.pi/(N_L*mpi))
    plt.ylabel("$\sqrt{s}$/$m_\pi$")
    plt.xlabel("$N_L$")
    plt.grid()
    xarr = np.linspace(min(NL_arr)-1, max(NL_arr)+1)
    pipi1arr = [sqrt_s_pi_L(x,1) for x in xarr]
    print(pipi1arr)
    pipi2arr = [sqrt_s_pi_L(x,2) for x in xarr]
    plt.axhline(2, color = "green", label = "pipi")
    plt.axhline(mrho/mpi, color = "orange", label = "rho")
    plt.plot(xarr, pipi1arr, ls = "dashdot", color = "green")
    plt.plot(xarr, pipi2arr, ls = "dashdot", color = "green")

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

    for i in range(len(sqrt_s_arr)):
        plt.errorbar(NL_arr[i], sqrt_s_arr[i], yerr = [[abs(sqrt_s_m_arr[i]-sqrt_s_arr[i]),],[abs(sqrt_s_p_arr[i]-sqrt_s_arr[i]),]],marker = marker(P_arr[i]), color = colorf(en_lv_arr[i]), capsize=3)

    plt.title("red exited, 001-o, 110-x")

    plt.legend()
    plt.savefig("sqrt_s_l.pdf")
    plt.show()

def plot_e():
    def en_L(N_L,p2):
        return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
    def en_rho_L(N_L,p2):
        return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi
    plt.ylabel("E/m$\pi$")
    plt.xlabel("N_L")
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
        plt.errorbar(NL_arr[i], en_arr[i], yerr = [[en_m_arr[i],],[en_p_arr[i],]],marker = marker(P_arr[i]), color = colorf(en_lv_arr[i]), capsize=3)

    plt.title("red exited, 001-o, 110-x")

    plt.legend()
    plt.savefig("en_l.pdf")
    plt.show()
    plt.clf()

def plot_pcot_PS():
    def en_L(N_L,p2):
        return 1+np.sqrt(1+p2*(2*np.pi/(N_L*mpi))**2)
    def en_rho_L(N_L,p2):
        return np.sqrt(mrho**2+p2*(2*np.pi/N_L)**2)/mpi
    plt.ylabel("E/m$\pi$")
    plt.xlabel("N_L")
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
        plt.errorbar(NL_arr[i], en_arr[i], yerr = [[en_m_arr[i],],[en_p_arr[i],]],marker = marker(P_arr[i]), color = colorf(en_lv_arr[i]), capsize=3)

    plt.title("red exited, 001-o, 110-x")

    plt.legend()
    plt.savefig("en_l.pdf")
    plt.show()
    plt.clf()


if __name__ == "__main__":
    plot_e()    
    # plot_sqrt_s()