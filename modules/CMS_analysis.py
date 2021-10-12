import numpy as np
from copy import deepcopy
from event_parser import data_set
import basic_plotter as plotter
import vmath
import physics as phys
import time

#Note: OS is implicit here
def list_choices():
    print("Set filter by inputting an array of choices")
    print("Event:0:trig")
    print("Event:1:lxy_cut")
    print("Pair:2:iso_cut")
    print("Pair:3:dphi_muonsDV_cut")
    print("Pair:4:dphi_muons_cut")
    print("Pair:5:pileup_cut")
    print("Pair:6:dxy_cut")

def print_cut_selection(arr):
    print("Cut choices")
    if(arr[0]):
        print("True:trig")
    else:
        print("False:trig")
    
    if(arr[1]):
        print("True:lxy_cut")
    else:
        print("False:lxy_cut")
    
    if(arr[2]):
        print("True:iso_cut")
    else:
        print("False:iso_cut")

    if(arr[3]):
        print("True:dphi_muonsDV_cut")
    else:
        print("False:dphi_muonsDV_cut")

    if(arr[4]):
        print("True:dphi_muons_cut")
    else:
        print("False:dphi_muons_cut")
    
    if(arr[5]):
        print("True:pileup_cut")
    else:
        print("False:pileup_cut")
    
    if(arr[6]):
        print("True:dxy_cut")
    else:
        print("False:dxy_cut")

def trigA(mu, amu):
    #Calculate muon deltaR
    muon_deltaR = phys.deltaR(mu, amu)
    cut_pT = (mu.pT > 4.) & (amu.pT > 4.)
    cut_deltaR = muon_deltaR < 1.2
    cut = cut_pT & cut_deltaR
    return cut

def trigB(mu, amu):
    muon_deltaR = phys.deltaR(mu, amu)
    cut_eta = (np.abs(mu.eta) < 1.4) & (np.abs(amu.eta) < 1.4)
    cut_deltaR = muon_deltaR < 1.4
    cut = cut_eta & cut_deltaR
    return cut

def trigC(mu, amu):
    c1 = mu.pT > 15.
    c2 = mu.pT > 7.
    c3 = amu.pT > 15.
    c4 = amu.pT > 7.

    cut = (c1 & c4) | (c2 & c3)  
    return cut

def trig(mu, amu):
    cA = trigA(mu, amu)
    cB = trigB(mu, amu)
    cC = trigC(mu, amu)
    cD = (mu.pT > 3.) & (amu.pT > 3.) & (np.abs(mu.eta) < 2.4) & (np.abs(amu.eta) < 2.4)

    cut = (cA | cB | cC) & cD
    return cut

def lxy_cut(mu):
    #Note: Both mu and amu have the same lxy
    lxy = phys.lxy(mu)
    cut = lxy < 110. #mm here, CMS paper is in cm
    return cut

def dphi_muonsDV_cut(mu, amu, n):
    max_dphi = 0.

    if(n == 1):
        max_dphi = 0.02
    elif(n == 2):
        max_dphi = 0.01
    else:
        print("Not a valid muon n")
    
    DV = deepcopy(mu.xi)
    DV[:,3] = 0.
    pmm = (mu.p + amu.p)
    pmm[:,3] = 0.
    dphi = vmath.vec4_dphi(DV, pmm)
    cut = dphi < max_dphi 
    return cut

def iso_cut(mu, amu, tr, jets, n): 
    R = 0.3
    mu_iso = phys.iso(mu, tr, R, 'track')
    amu_iso = phys.iso(amu, tr, R, 'track')
    if(n == 1):
        iso_cut = (mu_iso > 0.1*mu.pT) | (amu_iso > 0.1*amu.pT)
    elif(n == 2):
        iso_cut = (mu_iso > 0.2*mu.pT) | (amu_iso > 0.2*amu.pT)

    mu_jets = phys.iso(mu, jets, R, 'jet').astype(bool)
    amu_jets = phys.iso(amu, jets, R, 'jet').astype(bool)
    jet_cut = ~mu_jets | ~amu_jets

    cut = jet_cut | iso_cut
    return cut

def dphi_muons_cut(mu, amu):
    dphi = vmath.dphi(mu.phi, amu.phi)
    cut = dphi < 2.8
    return cut

def pileup_cut(mu, amu):
    pileup = phys.pileup(mu,amu)
    cut = pileup < 1.25
    return cut

def dxy_cut(sigma, mu, amu, n):
    ##See the bottom of CMS paper page 5
    sig = sigma # 1mm sigma
    mu_dxy = phys.dxy(mu)
    mu_lxy = phys.lxy(mu)
    amu_dxy = phys.dxy(amu)
    amu_lxy = phys.lxy(amu)
    
    min_vbl = 0.
    min_dxysig = 0.

    if(n == 1):
        min_vbl = 0.1
        min_dxysig = 2.
    elif(n == 2):
        min_vbl = 0.05
        min_dxysig = 1.
    else:
        print("Not a valid n")

    cut1 = (np.abs(mu_dxy)/sig) > min_dxysig
    
    m_12 = vmath.vec4_invar(mu.p + amu.p) #Check
    p_sum = deepcopy(mu.p) + deepcopy(amu.p)
    p_sum[:,3] = 0.
    pT_12 = vmath.vec4_norm(p_sum)
    mu_vbl = mu_dxy/(mu_lxy*m_12/pT_12)
    amu_vbl = amu_dxy/(amu_lxy*m_12/pT_12)
    cut2 = (mu_vbl > min_vbl) & (amu_vbl > min_vbl)

    cut = cut1 & cut2
    return cut

def event_cut(arr, mu, amu):
    cut = np.ones(mu.size()).astype(bool)
    
    if (arr[0]):
        cut = cut & trig(mu, amu) 
    if (arr[1]):
        cut = cut & lxy_cut(mu)
    return cut

def pair_cut(sigma, arr, mu, amu, tracks, jets, n):
    cut = np.ones(mu.size()).astype(bool)

    if (arr[2]):
        cut = cut & iso_cut(mu, amu, tracks, jets, n)
    if (arr[3]):
        cut = cut & dphi_muonsDV_cut(mu, amu, n)
    if (arr[4]):
        cut = cut & dphi_muons_cut(mu, amu)
    if (arr[5]):
        cut = cut & pileup_cut(mu, amu)
    if (arr[6]):
        cut = cut & dxy_cut(sigma, mu, amu, n)

    return cut

def isolate_repeats(data):#Should move this function into evt_parser
    #Find all pairs of repeats
    repeats = np.zeros(dp.size(),dtype=bool)

    for i in range(data.size()-1):
        if data.event_index[i] == data.event_index[i+1]:
            repeats[i] = True
            repeats[i+1] = True
            i = i+2

    return repeats

def prep_dataset(my_set):
    #Note that this function will need to be edited again if you allow decays to non-muons

    #Split out the dark photons
    dark_photons = my_set.get_pid(999999)

    #Get just muons and antimuons, + the set of aligned dark photons
    daughters = my_set.get_type(1)
    muons = daughters.get_pid(13)
    antimuons = daughters.get_pid(-13)

    #Get other entries 
    tracks = my_set.get_type(2)
    jets = my_set.get_type(3)
    jets.eta = phys.calc_pseudo_rap(jets)

    return dark_photons, muons, antimuons, tracks, jets

def get_cms_eff(sigma, dark_photons, muons, antimuons, tracks, jets, arr):
    #Function splits out full cut and re-ordered darkphotons, muons, antimuons
    #For efficiency correlogram

    #Load dataset from filename    
    start_time = time.time()
    
    init = muons.size()
    n = 1 #First/second muon pair - for non-repeat events, all are first
    
    #Find all repeat event entries and split sets accordingly
    #r_entries = dark_photons.find_repeats()
    r_entries = muons.find_repeats()
    r_dp, s_dp = dark_photons.split_set(r_entries)
    r_mu, s_mu = muons.split_set(r_entries)
    print("r: " + str(r_mu.size()))
    print("s: " + str(s_mu.size()))
    r_amu, s_amu = antimuons.split_set(r_entries)

    #Run all single entries through cuts 
    s_evt_cut = event_cut(arr, s_mu, s_amu)
    s_pair_cut = pair_cut(sigma, arr, s_mu, s_amu, tracks, jets, 1)
    s_cut = s_evt_cut & s_pair_cut
    
    final_dp, final_mu, final_amu = s_dp, s_mu, s_amu
    final_cut = s_cut
    final_evt_cut = s_evt_cut
    
    #Split double events along even and odd entries
    eo_split = r_dp.create_even_odd_split()
    er_dp, or_dp = r_dp.split_set(eo_split)
    er_mu, or_mu = r_mu.split_set(eo_split)
    er_amu, or_amu = r_amu.split_set(eo_split)
    
    #Create cut where one of even or odd passes the event cut
    #The other is included even if it doesn't pass cut criteria
    er_cut = event_cut(arr, er_mu, er_amu)
    or_cut = event_cut(arr, or_mu, or_amu)
    evt_cut = er_cut | or_cut
    
    #Split muon/antimuon set into closer/farther from beamline
    a = vmath.vec3_norm(er_mu.xi)
    b = vmath.vec3_norm(or_mu.xi)
    mask = a > b
    ar_mu = er_mu.get(mask).append(or_mu.get(~mask))
    br_mu = er_mu.get(~mask).append(or_mu.get(mask))
    ar_amu = er_amu.get(mask).append(or_amu.get(~mask))
    br_amu = er_amu.get(~mask).append(or_amu.get(mask))
    ar_dp = er_dp.get(mask).append(or_dp.get(~mask))
    br_dp = er_dp.get(~mask).append(or_dp.get(mask))
    evt_cutA = np.append(evt_cut[mask], evt_cut[~mask])
    evt_cutB = np.append(evt_cut[~mask], evt_cut[mask])

    #Next run er events through with stringent conditions, 
    #and or through with loosened conditions
    ar_pair_cut = pair_cut(sigma, arr, ar_mu, ar_amu, tracks, jets, 1)
    br_pair_cut = pair_cut(sigma, arr, br_mu, br_amu, tracks, jets, 2)
    final_cut = np.append(final_cut, (ar_pair_cut & evt_cutA))
    final_cut = np.append(final_cut, (br_pair_cut & evt_cutB))
    final_evt_cut = np.append(final_evt_cut, evt_cutA)
    final_evt_cut = np.append(final_evt_cut, evt_cutB)
    
    #Finally, append all particles
    final_dp = final_dp.append(ar_dp)
    final_dp = final_dp.append(br_dp)
    final_mu = final_mu.append(ar_mu)
    final_mu = final_mu.append(br_mu)
    final_amu = final_amu.append(ar_amu)
    final_amu = final_amu.append(br_amu)
    
    print("Final: " + str(final_cut[final_cut].shape[0]) + " muon pairs")
    final = final_cut[final_cut].shape[0]
    print("Ratio: " + str(final/init))
    ratio = final/init
    print("Analyze time: --- %s seconds ---" % (time.time()-start_time))

    return final_dp, final_mu, final_amu, final_cut, final_evt_cut
