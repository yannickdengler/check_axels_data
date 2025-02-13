import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
import h5py



# data = np.transpose(np.genfromtxt("data/axel_data.dat"))
data = np.transpose(np.genfromtxt("data/corrfitter.dat"))

print(data.shape)

mpi = 0.38649
mrho = 0.5494



# should work with m1=m2=1 as it only is needed for mass degenerate case



NL_arr = data[0]
P_arr = data[4]
Pvec_arr = np.transpose(data[1:4])
P_m_arr = [np.sqrt(np.dot(Pvec_arr[i],Pvec_arr[i]))*2*np.pi/(mpi*NL_arr[i]) for i in range(len(Pvec_arr))]
en_lv_arr = data[5]
en_arr = data[6]/mpi
en_p_arr = data[7]/mpi
en_m_arr = data[8]/mpi

def save_to_hdf(res,res_sample, filename):
    with h5py.File("output/hdf5/"+filename+".hdf5","w") as hfile:
        for key, val in res.items():
            hfile.create_dataset("orig_"+key, data = val)
        for key, val in res_sample.items():
            hfile.create_dataset("sample_"+key, data = val)

def read_from_hdf(filename):
    res, res_tmp = [{},{}]
    with h5py.File("output/hdf5/"+filename+".hdf5","r") as hfile:
        for key in hfile.keys():
            if key[:4] == "orig":
                res[key[5:]] = hfile[key][()]
            if key[:4] == "samp":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def result_sampled(beta,m0,N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,en_lv,mpi,num_gaussian=50):
    res = get_rizz(E_pipi,N_L,dvec,en_lv,mpi)
    res_sample = {}
    for key in res.keys():
        res_sample[key] = []

    for i in tqdm(range(num_gaussian)):
        if i < num_gaussian//2:
            E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
        else:
            E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
    # for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_gaussian)):
        res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,en_lv,mpi)
        for key, val in res_tmp.items():
            res_sample[key].append(val)
    res["beta"],res["m_1"],res["m_2"],res_sample["beta"],res_sample["m_1"],res_sample["m_2"] = [beta,m0,m0,beta,m0,m0]
    return res, res_sample

def Ecm(E, P):
    return np.sqrt(E**2-P**2)

def pstar(Ecm):
    return np.sqrt(Ecm**2/4-1)

def cot_delta_001(q2,N_L):
    return wlm(0,0,0,0,1,1,1,q2,int(N_L)) + 2*wlm(2,0,0,0,1,1,1,q2,int(N_L))
def cot_delta_110(q2,N_L):
    return wlm(0,0,1,1,0,1,1,q2,int(N_L)) - wlm(2,0,1,1,0,1,1,q2,int(N_L)) + np.sqrt(3/2) * (wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L)))

def cot_delta_mom(dvec):
    if list(dvec) == [0,0,1]:
        return cot_delta_001
    elif list(dvec) == [1,1,0]:
        return cot_delta_110
    else:
        print("wrong momentum")
        exit()

def get_rizz(E_pipis, N_Ls, dvecs, en_lvs, mpi):
    result = {}

    E_prime = []
    d = []
    L_prime = []
    P_prime = []
    E_cm_prime = []
    s_prime = []
    pstar_prime = []
    pstar_2_prime = []
    q = []
    q2 = []
    cot_PS = []
    tan_PS = []
    P_cot_PS_prime = []
    P3_cot_PS_prime = []

    for i in range(len(E_pipis)):
        E_prime_tmp = E_pipis[i]/mpi
        d_tmp = np.sqrt(np.dot(dvecs[i],dvecs[i]))
        L_prime_tmp = N_Ls[i]*mpi                               # L*mpi (removes lattice const)
        P_prime_tmp = 2*np.pi*d_tmp/L_prime_tmp
        
        if E_prime_tmp**2 - P_prime_tmp**2 < 4:
            E_cm_prime_tmp,s_prime_tmp,pstar_prime_tmp,pstar_2_prime_tmp,q_tmp,q2_tmp,cot_PS_tmp,tan_PS_tmp,P_cot_PS_prime_tmp,P3_cot_PS_prime_tmp =  [0,0,0,0,0,0,0,0,0,0]
        else:
            E_cm_prime_tmp = Ecm(E_prime_tmp,P_prime_tmp)
            s_prime_tmp = E_cm_prime_tmp**2
            pstar_prime_tmp = pstar(E_cm_prime_tmp)
            pstar_2_prime_tmp = pstar_prime_tmp**2
            q_tmp = pstar_prime_tmp*L_prime_tmp/(2*np.pi)
            q2_tmp = q_tmp**2
            cot_PS_tmp = cot_delta_mom(dvecs[i])(q2_tmp,N_Ls[i])
            tan_PS_tmp = 1/cot_PS_tmp
            P_cot_PS_prime_tmp = cot_PS_tmp*pstar_prime_tmp
            P3_cot_PS_prime_tmp = cot_PS_tmp*pstar_prime_tmp**3
        
        E_prime.append(E_prime_tmp)
        d.append(d_tmp)
        L_prime.append(L_prime_tmp)
        E_cm_prime.append(E_cm_prime_tmp)
        s_prime.append(s_prime_tmp)
        pstar_prime.append(pstar_prime_tmp)
        pstar_2_prime.append(pstar_2_prime_tmp)
        q.append(q_tmp)
        q2.append(q2_tmp)
        cot_PS.append(cot_PS_tmp)
        tan_PS.append(tan_PS_tmp)
        P_cot_PS_prime.append(P_cot_PS_prime_tmp)
        P3_cot_PS_prime.append(P3_cot_PS_prime_tmp)

    result["level"] = en_lvs
    result["E_pure"] = E_pipis
    result["E_prime"] = E_prime
    result["dvecs"] = dvecs
    result["d"] = d
    result["N_Ls"] = N_Ls
    result["E_cm_prime"] = E_cm_prime
    result["s_prime"] = s_prime
    result["pstar_prime"] = pstar_prime
    result["pstar_2_prime"] = pstar_2_prime
    result["L_prime"] = L_prime
    result["q"] = q
    result["q2"] = q2
    result["cot_PS"] = cot_PS
    result["tan_PS"] = tan_PS
    result["PS"] = np.arctan(tan_PS)
    result["P_cot_PS_prime"] = P_cot_PS_prime
    result["P3_cot_PS_prime"] = P3_cot_PS_prime
    return result


def calc_PS():
    NL_arr = data[0]
    dvec_arr = np.transpose(data[1:4])
    en_lv_arr = data[5]
    en_arr = data[6]
    en_m_arr = data[7]
    en_p_arr = data[8]
    mpi = 0.38649

    res, res_sampled = result_sampled(6.9, -0.92, NL_arr, en_arr, en_m_arr, en_p_arr, dvec_arr, en_lv_arr, mpi)
    save_to_hdf(res, res_sampled, "PS_69_092")

if __name__ == "__main__":
    calc_PS()