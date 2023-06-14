import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

def basin_entropy(w):
    # use masked array (ma) to avoid problems with zero
    return -np.dot(w,np.ma.log2(w).filled(0.))

def entropy_to_relative_basin_sizes(h_tilde,r,n,tol=1e-3):
    
    h_tilde_min = (r-1)*n/2**n -(1-(r-1)/2**n)*np.log2((1-(r-1)/2**n))
    
    if h_tilde > np.log2(r): 
        raise ValueError("Basin entropy cannot be larger than log2(r)")
        
    if h_tilde < h_tilde_min:
        raise Exception("h_tilde below h_tilde_min = "+str(h_tilde_min))
        
    func = lambda w: (basin_entropy(w) - h_tilde)**2
    w0 = [ (i+1)/(r*(r+1)/2) for i in range(r) ]
    normConstraint = scipy.optimize.LinearConstraint(np.ones(r),lb=1.,ub=1.)
    bound = scipy.optimize.Bounds(lb=1./2**n,ub=1.)
    solution = scipy.optimize.minimize(func,w0,constraints=normConstraint,bounds=bound)
    if solution.success and np.sqrt(solution.fun) < tol:
        return solution.x
    else:
        raise Exception("Solution not found within given tolerance")

def entropy_to_basin_sizes(h_tilde,r,n):
	"""
	To do: Add checks to this!!!
	"""
    w = entropy_to_relative_basin_sizes(h_tilde,r,n)
    w_int = [round(w[_]*2**n) for _ in range(len(w))]
    w_int[-1] = w_int[-1] +  2**n - np.sum(w_int)
    return np.array(w_int)



