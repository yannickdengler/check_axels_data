import h5py
import numpy as np

axel_data = np.transpose(np.genfromtxt("data/axel_data.dat"))
NL_arr = [int(x) for x in axel_data[0]]

print(axel_data.shape)

en_lv = []
en_lv_ep = []
en_lv_em = []
P_arr = []

for i in range(len(axel_data)):
    Ptmp = str(axel_data[1:4,i])
    print(Ptmp)
    print(P_arr)
    if not(Ptmp in P_arr):
        P_arr.append(Ptmp)
        en_lv.append([])
        en_lv_ep.append([])
        en_lv_em.append([])
    en_lv[P_arr.index(Ptmp)].append(axel_data[6,i])
    en_lv_ep[P_arr.index(Ptmp)].append(axel_data[7,i])
    en_lv_em[P_arr.index(Ptmp)].append(axel_data[8,i])

print(en_lv)



with h5py.open("h5py_fit_results.hdf5", "w") as f:
    f.create_dataset("isospin_channel", data="1")
    f.create_dataset("gauge_group", data = "Sp(4)")
    f.create_dataset("beta", data = 6.9)
    f.create_dataset("m_1", data = -0.92)
    f.create_dataset("m_2", data = -0.92)
    f.create_dataset("N_Ls", data = NL_arr)
    f.create_dataset("correlators", data = Correlators)