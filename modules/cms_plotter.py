import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from vhepmc_evt_parser import data_set
import vmath as vm
import physics as phys
from math import *


##First portion is about efficiency plots
##Functions
def get_eff_divs(edges, arr):
    #Initialize an empty array
    arr_divs = np.array([])

    #Append to arr_divs values of arr for which edges[i] < arr and edges[i+1] > arr
    for i in range(edges.shape[0]-1):
        arr_divs = np.append(arr_divs,(edges[i] < arr) & (arr < edges[i+1]), axis=0)

    #Make sure the array is full of booleans
    arr_divs = np.array(arr_divs).astype(bool)

    #Reshape such that function returns an array of arrays matching each edges[i] < arr & etc. cond.
    arr_divs = arr_divs.reshape( (edges.shape[0]-1), int(arr_divs.shape[0]/(edges.shape[0]-1)))
    return arr_divs

def get_eff(arr1_divs, arr2_divs, cms_cut):
   
    #Initialize an empty array -- note that  
    eff = np.zeros((arr1_divs.shape[0], arr2_divs.shape[0])) 
    
    #Going across rows
    for i in range(arr1_divs.shape[0]):
        #Going across columns
        for j in range(arr2_divs.shape[0]):
            #Entries that match both conditions are true
            split = (arr1_divs[i] & arr2_divs[j])
            if(split[split].shape[0] == 0):
                #If no entries meet this category, fill with nan
                eff[i][j] = np.nan
            else:
                #Efficiency = (number that passed cms_cut+split)/(number that fit in split)
                rem = cms_cut & split
                eff[i][j] = rem[rem].shape[0]/split[split].shape[0]

    eff = np.flip(eff, axis=0)
    
    return eff

def get_n_eff(arr1_divs, arr2_divs, cms_cut):
   
    #Initialize an empty array with as many entries along 
    eff = np.zeros((arr1_divs.shape[0], arr2_divs.shape[0])) 
    tots = np.zeros((arr1_divs.shape[0], arr2_divs.shape[0])) 
    #Going across rows
    for i in range(arr1_divs.shape[0]):
        #Going across columns
        for j in range(arr2_divs.shape[0]):
            #Entries that match both conditions are true
            split = (arr1_divs[i] & arr2_divs[j])
            if(split[split].shape[0] == 0):
                #If no entries meet this category, fill with nan
                eff[i][j] = np.nan
                tots[i][j] = 0.
            else:
                #Efficiency = (number that passed cms_cut+split)/(number that fit in split)
                rem = cms_cut & split
                eff[i][j] = rem[rem].shape[0]
                tots[i][j] = split[split].shape[0]

    eff = np.flip(eff, axis=0)
    tots = np.flip(tots, axis=0)

    return eff, tots

def create_eff_plot(edges1, edges2, arr1, arr2, cut, save=False, name="my_name.png"):
    #Index
    ind = []
    index_vbl = "l_xy"
    
    for i in range(edges1.shape[0]-1):
        ind.append(str(edges1[i]) + " < " + index_vbl + " < " + str(edges1[i+1]))
    
    ind.reverse()

    #Columns
    col_vbl = 'pT'
    col = []
    for i in range(edges2.shape[0]-1):
        col.append(str(edges2[i]) + " < " + col_vbl + " < " + str(edges2[i+1]))


    #Efficiency
    eff = get_eff(arr1, arr2, cut)
    df = pd.DataFrame(eff, index=ind, columns=col)

    #Plot
    fig = plt.figure(figsize=(12,9))
    sns.heatmap(df,cmap='Blues', annot=True, vmin=0)
    plt.show()
    
    if(save):
        fig.savefig("fig_v1.png")

def create_eff_plot_v2(edges1, edges2, eff, plt_type, save=False, name="my_name.png"):
    #Index
    ind = []
    index_vbl = "l_xy"
    
    for i in range(edges1.shape[0]-1):
        ind.append(str(edges1[i]) + " < " + index_vbl + " < " + str(edges1[i+1]))
    
    ind.reverse()

    #Columns
    col_vbl = 'pT'
    col = []
    for i in range(edges2.shape[0]-1):
        col.append(str(edges2[i]) + " < " + col_vbl + " < " + str(edges2[i+1]))

    df = pd.DataFrame(eff, index=ind, columns=col)

    #Plot
    fig = plt.figure(figsize=(12,9))
    if(plt_type == 1):
        sns.heatmap(df,cmap='Blues', annot=True, vmin=0)
    if(plt_type == 2):
        sns.heatmap(df,cmap='YlOrBr', annot=True, vmin=0)
    plt.show()
    
    if(save):
        fig.savefig(name)
