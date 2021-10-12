import numpy as np
from math import *
from copy import deepcopy 

def vec4_dotprod(v1, v2):
    dot = v1[:,1]*v2[:,1] + v1[:,2]*v2[:,2] + v1[:,3]*v2[:,3]
    return dot

def vec4_norm(v):
    norm = np.sqrt(v[:,1]*v[:,1] + v[:,2]*v[:,2] + v[:,3]*v[:,3])
    return norm

def vec4_dang(v1, v2):
    ang = np.arccos(vec4_dotprod(v1,v2)/(vec4_norm(v1)*vec4_norm(v2)))
    return ang

def vec4_dphi(v1, v2):
    c_v1 = deepcopy(v1)
    c_v2 = deepcopy(v2)
    c_v1[:,3] = 0. #Note: operations like this need deepcopy b/c pass by ref.
    c_v2[:,3] = 0.

    ang = vec4_dotprod(c_v1,c_v2)/(vec4_norm(c_v1)*vec4_norm(c_v2))
    #Sometimes greater than or less than 1 - shouldn't be possible.
    #it's minimal, so I'm assuming its a rounding area.
    mask = (ang < 1.1) & (ang > 1.) 
    ang[mask] = 1.
    mask = (ang > -1.1) & (ang < -1.)
    ang[mask] = -1.
    
    ang = np.arccos(ang)
    return ang

def dphi(phi1, phi2):
    dphi = np.abs(phi1-phi2)
    mask = dphi > np.pi 
    dphi[mask] = 2*np.pi-dphi[mask] 
    return dphi

def vec4_err(v1, v2):
    v1[:,3] = 0.
    v2[:,3] = 0.
    ang = vec4_dotprod(v1,v2)/(vec4_norm(v1)*vec4_norm(v2))
    check = ang > 1.
    print(check[check].shape)
    return ang

def vec4_invar(v):
    return np.sqrt(v[:,0]*v[:,0] - v[:,1]*v[:,1] - v[:,2]*v[:,2] - v[:,3]*v[:,3])

def vec3_norm(v):
    norm = np.sqrt(v[:,0]*v[:,0] + v[:,1]*v[:,1] + v[:,2]*v[:,2])
    return norm

def vec3_dotprod(v1, v2):
    dot = v1[:,0]*v2[:,0] + v1[:,1]*v2[:,1] + v1[:,2]*v2[:,2]
    return dot
