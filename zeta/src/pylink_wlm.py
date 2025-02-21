import ctypes as ct
import sys

class c_double_complex(ct.Structure): 
    """complex is a c structure
    https://docs.python.org/3/library/ctypes.html#module-ctypes suggests
    to use ctypes.Structure to pass structures (and, therefore, complex)
    """
    _fields_ = [("real", ct.c_double),("imag", ct.c_double)]
    @property
    def value(self):
        return self.real+1j*self.imag # fields declared above

C_wlm = ct.CDLL(sys.path[0]+"/../out/get_wlm.so")

wlm = C_wlm.wlm_PL

wlm.argtypes = [ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_int]
wlm.restype = c_double_complex

def wlm_func(l,m,d1,d2,d3,m1,m2,q2,N_L,size=10):            # m1, m2 in lattice units (m*a)
    res = wlm(l,m,d1,d2,d3,m1,m2,q2,N_L,size).value
    return res.real, res.imag

def wlm_func_c(l,m,d1,d2,d3,m1,m2,q2,N_L,size=10):            # m1, m2 in lattice units (m*a)
    res = wlm(l,m,d1,d2,d3,m1,m2,q2,N_L,size).value
    return res

if __name__ == "__main__":
    m1=float(sys.argv[1])
    m2=float(sys.argv[2])
    q2=float(sys.argv[3])
    l=int(sys.argv[4])
    m=int(sys.argv[5])
    N_L=int(sys.argv[6])
    sz=int(sys.argv[7])
    d1=int(sys.argv[8])
    d2=int(sys.argv[9])
    d3=int(sys.argv[10])
    val = wlm_func(l,m,d1,d2,d3,m1,m2,q2,N_L,sz)
    print("%18.14f,%18.14f"%(val[0], val[1]))
