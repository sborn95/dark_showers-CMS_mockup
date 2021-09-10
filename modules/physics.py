import numpy as np
import vmath as vm
from math import *
from copy import deepcopy
from vhepmc_evt_parser import data_set

def comb_pT(d1, d2):
    p = d1.p + d2.p
    p[:,3] = 0
    pT = vm.vec4_norm(p)
    return pT

def deltaR(d1, d2):
    d_eta = d1.eta-d2.eta
    d_phi = vm.vec4_dphi(d1.p, d2.p)
    return np.sqrt(d_eta*d_eta + d_phi*d_phi)

def lxy(d):
    vxy = deepcopy(d.xi)
    vxy[:,3] = 0.
    lxy = vm.vec4_norm(vxy)
    return lxy

def dxy(d):
    a = deepcopy(d.xi)[:,1:4] #Point on the line
    a[:,2] = 0. #Just looking at x-y plane projection
    p = deepcopy(d.p)[:,1:4]
    p[:,2] = 0.
    normp = np.zeros(d.size()*3).reshape(d.size(),3)
    temp = vm.vec3_norm(p)
    normp[:,0] = temp
    normp[:,1] = temp
    normp[:,2] = temp
    n = p/(normp)
    X = deepcopy(d.xPV)[:,1:4]
    proj = np.zeros(d.size()*3).reshape(d.size(),3)
    temp = vm.vec3_dotprod((X-a),n)
    proj[:,0] = temp
    proj[:,1] = temp
    proj[:,2] = temp
    dxy = vm.vec3_norm((X-a) - (proj*n))
    return dxy

def pileup(d1, d2):
    deltaPhi = vm.vec4_dphi(d1.p,d2.p)
    deltaEta = np.abs(d1.eta-d2.eta)
    pileup = np.log10(deltaEta/deltaPhi)
    return pileup

def calc_pseudo_rap(d):
    cosA = d.p[:,3]/vm.vec4_norm(d.p)
    A = np.arccos(cosA)
    pr = -np.log(np.tan(A/2))
    return pr

def iso(mu, jtr, R, typ):
    mu_i = data_set()
    mu_iso = np.array([])
    found_jet = np.array([])
    
    for i in range(mu.size()):
        #Get all tracks associated with the event index of muon i
        mu_i.clear_data()
        mask = jtr.event_index == mu.event_index[i]
        n = mask[mask].shape[0]
        jtr_i = jtr.get(mask)

        mu_i = mu.copy_single(4, n, i)
        mu_i.eta = mu.copy_single(6, n, i).eta

        #Computer deltaR
        dR = deltaR(jtr_i, mu_i)
        found = dR < R

        if(typ == 'jet'):
            if (found[found].shape[0] > 0): found_jet = np.append(found_jet, True)
            else: found_jet = np.append(found_jet, False)

        elif (typ == 'track'):
            cone_jtr = jtr_i.get(dR <= R)
            iso = np.sum(cone_jtr.pT)
            mu_iso = np.append(mu_iso, iso)

    if(typ == 'jet'): return found_jet
    elif (typ == 'track'): return mu_iso
    else:
        print("Not a valid typ entry")
        return 0


