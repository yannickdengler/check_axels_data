import h5py
import numpy as np

zz = complex(123, 3223)
zznp = np.complex(123, 3223)


print(type(zz))
# print(type(zznp))


with h5py.File("test_complex.hdf5", "w") as f:
    f.create_dataset("complex_numer", data=zz)

with h5py.File("test_complex.hdf5", "r") as f:
    print(f["complex_numer"][()])
    print(type(f["complex_numer"][()]))