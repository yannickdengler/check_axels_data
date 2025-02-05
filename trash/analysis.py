import numpy as np




data = np.genfromtxt("data/axel_data.csv")

print(data.shape)

data_b69m092 = data[:5]

print(data_b69m092.shape)

energy_levels_b69_m092 = data_b69m092[::,6:]


print(energy_levels_b69_m092.shape)

print(energy_levels_b69_m092)

NL_arr = [16,16,14,24,24]
# P_arr = [[0,0,1],[1,1,0],[0,0,1],[0,0,1],[1,1,0]]
P_arr = [1,2,1,1,2]


