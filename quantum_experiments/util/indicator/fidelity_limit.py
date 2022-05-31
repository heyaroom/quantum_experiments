import numpy as np

def single_fidelity(gt,t1,t2):
    f = (0.5+1/6.*np.exp(-gt/t1)+1/3.*np.exp(-gt/t2))
    return f
    
def two_fidelity(f1,f2):
    f = 0.2+0.05*(6*f1-2)*(6*f2-2)
    return f