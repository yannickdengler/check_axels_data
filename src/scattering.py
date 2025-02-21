import numpy as np
import matplotlib.pyplot as plt
from pylink_wlm import wlm_func_c as wlm
from tqdm import tqdm
import h5py


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
            if key[:5] == "orig_":
                res[key[5:]] = hfile[key][()]
            if key[:7] == "sample_":
                res_tmp[key[7:]] = hfile[key][()]
    return res, res_tmp

def result_sampled(info,N_L,E_pipi,E_pipi_em,E_pipi_ep,dvec,mpi,num_resample=50, resampling = "gauss"):
    res = {}
    res_sample = {}
    res["resampling"] = resampling
    res["num_resample"] = num_resample
    res["mpi"] = mpi
    res["N_L"] = N_L
    res["L_prime"] = N_L*mpi
    res["dvec"] = dvec
    res["d"] = [np.sqrt(np.dot(dvec[i],dvec[i])) for i in range(len(dvec))]
    res["d2"] = [np.dot(dvec[i],dvec[i]) for i in range(len(dvec))]

    for key, val in info.items():
        res[key] = val
    for key, val in get_rizz(E_pipi,N_L,dvec,mpi).items():
        res[key] = val
    for key in res.keys():
        res_sample[key] = []

    if resampling == "gauss":
        for i in tqdm(range(num_resample)):
            if i < num_resample//2:
                E_pipi_tmp = E_pipi+abs(np.random.normal(0,E_pipi_ep))
            else:
                E_pipi_tmp = E_pipi-abs(np.random.normal(0,E_pipi_em))
    elif resampling == "lin":
        for E_pipi_tmp in tqdm(np.linspace(E_pipi-E_pipi_em,E_pipi+E_pipi_ep, num_resample)):
            res_tmp = get_rizz(E_pipi_tmp,N_L,dvec,mpi)
            for key, val in res_tmp.items():
                res_sample[key].append(val)
    # for key, val in res.items():
    #     print(key, val)
    
    # res["beta"],res["m_1"],res["m_2"], res["mpi"], res["mrho"], res["en_lv"],res["resampling"], res["num_resample"] = [beta,m0,m0,mpi,mrho, en_lv,resampling,num_resample]
    return res, res_sample


def cot_delta_000(q2,N_L):
    return wlm(0,0,0,0,0,1,1,q2,int(N_L))
def cot_delta_001(q2,N_L):
    return wlm(0,0,0,0,1,1,1,q2,int(N_L)) + 2*wlm(2,0,0,0,1,1,1,q2,int(N_L))
def cot_delta_110(q2,N_L):
    # print(complex(0,1) * (wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L))))
    # print(-2*np.imag(wlm(2,2,1,1,0,1,1,q2,int(N_L))))
    # third = wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L))
    # if np.absolute(third) > 1e-10:
    #     print(third, q2, N_L)
    return wlm(0,0,1,1,0,1,1,q2,int(N_L)) - wlm(2,0,1,1,0,1,1,q2,int(N_L)) + np.sqrt(3/2) * complex(0,1) * (wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L)))
    # return wlm(0,0,1,1,0,1,1,q2,int(N_L)) - wlm(2,0,1,1,0,1,1,q2,int(N_L)) - np.sqrt(3/2) * (wlm(2,2,1,1,0,1,1,q2,int(N_L)) - wlm(2,-2,1,1,0,1,1,q2,int(N_L)))

def cot_delta_mom(dvec):
    if list(dvec) == [0,0,0]:
        return cot_delta_000
    if list(dvec) == [0,0,1]:
        return cot_delta_001
    elif list(dvec) == [1,1,0]:
        return cot_delta_110
    else:
        print("wrong momentum")
        exit()

def Ecm_prime(E_prime, P_prime):
    return np.sqrt(E_prime**2-P_prime**2)

def pstar_prime(Ecm_prime):
    return np.sqrt(Ecm_prime**2/4-1)

def Ecm_lat_disp(E, Pvec):
    return np.arccosh(np.cosh(E)-2*(np.sin(Pvec[0]/2)**2+np.sin(Pvec[1]/2)**2+np.sin(Pvec[2]/2)**2))

def pstar_lat_disp(Ecm, mpi):
    return 2*np.arcsin(np.sqrt(0.5*(np.cosh(Ecm/2)-np.cosh(mpi))))

def get_rizz(E_pipis, N_Ls, dvecs, mpi):
    result = {}
    key_list = ["En","En_prime","E_cm","E_cm_prime","E_cm_ld","E_cm_ld_prime","s","s_prime","s_ld","s_ld_prime","pstar","pstar_prime","pstar_ld","pstar_ld_prime","p2star","p2star_prime","p2star_ld","p2star_ld_prime","q","q_ld","q2","q2_ld","cot_PS","cot_PS_ld","tan_PS","tan_PS_ld","PS","PS_ld", "p3cotPS", "p3cotPS_prime", "p3cotPS_ld", "p3cotPS_ld_prime", "p3cotPS_Ecm", "p3cotPS_Ecm_prime", "p3cotPS_Ecm_ld", "p3cotPS_Ecm_ld_prime"]

    for key in key_list:
        result[key] = []

    for i in range(len(E_pipis)):
        result["En"].append(E_pipis[i])
        Pvec = 2*np.pi*dvecs[i]/N_Ls[i]
        P_prime = 2*np.pi*np.sqrt(np.dot(dvecs[i],dvecs[i]))/(N_Ls[i]*mpi)
        En_prime = E_pipis[i]/mpi
        result["En_prime"].append(En_prime)
        if En_prime**2 - P_prime**2 < 4:
            for key in key_list[2:]:
                result[key].append(0)
        else:
            result["E_cm_prime"].append(Ecm_prime(E_pipis[i]/mpi,P_prime))
            result["E_cm"].append(result["E_cm_prime"][i]*mpi)
            result["s_prime"].append(result["E_cm_prime"][i]**2)
            result["s"].append(result["E_cm"][i]**2)
            result["E_cm_ld"].append(Ecm_lat_disp(E_pipis[i],Pvec))
            result["E_cm_ld_prime"].append(result["E_cm_ld"][i]/mpi)
            result["s_ld"].append(result["E_cm_ld"][i]**2)
            result["s_ld_prime"].append(result["E_cm_ld_prime"][i]**2)
            result["pstar_prime"].append(np.sqrt(result["s_prime"][i]/4-1))
            result["p2star_prime"].append(result["pstar_prime"][i]**2)
            result["pstar"].append(result["pstar_prime"][i]*mpi)
            result["p2star"].append(result["p2star_prime"][i]*mpi**2)
            result["pstar_ld"].append(pstar_lat_disp(result["E_cm_ld"][i],mpi))
            result["pstar_ld_prime"].append(result["pstar_ld"][i]/mpi)
            result["p2star_ld"].append(result["pstar_ld"][i]**2)
            result["p2star_ld_prime"].append(result["pstar_ld_prime"][i]**2)
            q2 = result["pstar"][i]*N_Ls[i]/(2*np.pi)
            q2_ld = result["pstar_ld"][i]*N_Ls[i]/(2*np.pi)
            result["q2"].append(q2)
            result["q2_ld"].append(q2_ld)
            cot_PS = cot_delta_mom(dvecs[i])(q2, N_Ls[i])
            cot_PS_ld = cot_delta_mom(dvecs[i])(q2_ld, N_Ls[i])
            result["cot_PS"].append(cot_PS)
            result["cot_PS_ld"].append(cot_PS_ld)
            result["tan_PS"].append(1/cot_PS)
            result["tan_PS_ld"].append(1/cot_PS_ld)
            PS = 360*np.arctan(1/cot_PS)/(2*np.pi)
            result["PS"].append(complex(PS.real%180,PS.imag%180))
            PS_ld = 360*np.arctan(1/cot_PS_ld)/(2*np.pi)
            result["PS_ld"].append(complex(PS_ld.real%180,PS_ld.imag%180))
            result["p3cotPS"].append(result["pstar"][i]**3)
            result["p3cotPS_prime"].append(result["pstar_prime"][i]**3)
            result["p3cotPS_ld"].append(result["pstar_ld"][i]**3)
            result["p3cotPS_ld_prime"].append(result["pstar_ld_prime"][i]**3)
            result["p3cotPS_Ecm"].append(result["pstar"][i]**3/result["E_cm"][i])
            result["p3cotPS_Ecm_prime"].append(result["pstar_prime"][i]**3/result["E_cm_prime"][i])
            result["p3cotPS_Ecm_ld"].append(result["pstar_ld"][i]**3/result["E_cm_ld"][i])
            result["p3cotPS_Ecm_prime"].append(result["pstar_ld_prime"][i]**3/result["E_cm_ld_prime"][i])
    return result


# def get_rizz(E_pipis, N_Ls, dvecs, mpi):
#     result = {}

#     E_prime = []
#     d = []
#     L_prime = []
#     P_prime = []
#     E_cm_prime = []
#     s_prime = []
#     pstar_prime = []
#     pstar_2_prime = []
#     q = []
#     q2 = []
#     cot_PS = []
#     tan_PS = []
#     PS = []
#     P_cot_PS_prime = []
#     P3_cot_PS_prime = []

#     for i in range(len(E_pipis)):
#         E_prime_tmp = E_pipis[i]/mpi
#         d_tmp = np.sqrt(np.dot(dvecs[i],dvecs[i]))
#         L_prime_tmp = N_Ls[i]*mpi                               # L*mpi (removes lattice const)
#         P_prime_tmp = 2*np.pi*d_tmp/L_prime_tmp
        
#         if E_prime_tmp**2 - P_prime_tmp**2 < 4:
#             E_cm_prime_tmp,s_prime_tmp,pstar_prime_tmp,pstar_2_prime_tmp,q_tmp,q2_tmp,cot_PS_tmp,tan_PS_tmp,PS_tmp,P_cot_PS_prime_tmp,P3_cot_PS_prime_tmp =  [0,0,0,0,0,0,0,0,0,0,0]
#         else:
#             E_cm_prime_tmp = Ecm(E_prime_tmp,P_prime_tmp)
#             s_prime_tmp = E_cm_prime_tmp**2
#             pstar_prime_tmp = pstar(E_cm_prime_tmp)
#             pstar_2_prime_tmp = pstar_prime_tmp**2
#             q_tmp = pstar_prime_tmp*L_prime_tmp/(2*np.pi)
#             q2_tmp = q_tmp**2
#             cot_PS_tmp = cot_delta_mom(dvecs[i])(q2_tmp,N_Ls[i])
#             tan_PS_tmp = 1/cot_PS_tmp
#             PS_tmp = complex((360*np.arctan(tan_PS_tmp)/(2*np.pi)).real%180,(360*np.arctan(tan_PS_tmp)/(2*np.pi)).imag%180)
#             P_cot_PS_prime_tmp = cot_PS_tmp*pstar_prime_tmp
#             P3_cot_PS_prime_tmp = cot_PS_tmp*pstar_prime_tmp**3
        
#         E_prime.append(E_prime_tmp)
#         d.append(d_tmp)
#         L_prime.append(L_prime_tmp)
#         E_cm_prime.append(E_cm_prime_tmp)
#         s_prime.append(s_prime_tmp)
#         pstar_prime.append(pstar_prime_tmp)
#         pstar_2_prime.append(pstar_2_prime_tmp)
#         q.append(q_tmp)
#         q2.append(q2_tmp)
#         cot_PS.append(cot_PS_tmp)
#         tan_PS.append(tan_PS_tmp)
#         PS.append(PS_tmp)
#         P_cot_PS_prime.append(P_cot_PS_prime_tmp)
#         P3_cot_PS_prime.append(P3_cot_PS_prime_tmp)

#     result["mpi"] = [mpi,]
#     result["level"] = en_lvs
#     result["E_pure"] = E_pipis
#     result["E_prime"] = E_prime
#     result["dvecs"] = dvecs
#     result["d"] = d
#     result["d2"] = d
#     result["N_Ls"] = N_Ls
#     result["E_cm_prime"] = E_cm_prime
#     result["s_prime"] = s_prime
#     result["pstar_prime"] = pstar_prime
#     result["pstar_2_prime"] = pstar_2_prime
#     result["L_prime"] = L_prime
#     result["q"] = q
#     result["q2"] = q2
#     result["cot_PS"] = cot_PS
#     result["tan_PS"] = tan_PS
#     result["PS"] = PS
#     result["P_cot_PS_prime"] = P_cot_PS_prime
#     result["P3_cot_PS_prime"] = P3_cot_PS_prime
#     return result

def calc_PS(name):
    info={}
    data = np.transpose(np.genfromtxt("data/%s.dat"%name))
    NL_arr = data[0]
    dvec_arr = np.transpose(data[1:4])
    en_lv_arr = data[4]
    en_arr = data[5]
    en_m_arr = data[6]
    en_p_arr = data[7]
    mpi = data[8][0]
    mrho = data[9][0]
    info["beta"],info["m_1"],info["m_2"], info["mrho"], info["en_lv"] = [6.9,-0.92,-0.92,mrho,en_lv_arr]

    res, res_sampled = result_sampled(info,NL_arr, en_arr, en_m_arr, en_p_arr, dvec_arr, mpi, resampling="lin")
    save_to_hdf(res, res_sampled, name)

if __name__ == "__main__":
    calc_PS("PS_69_092")
    # calc_PS("Plymouth")
    # calc_PS("Lang_Prelovsek")