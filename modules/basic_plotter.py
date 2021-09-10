import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from math import *

def plot_any(name, arr, nbins=50, rng=(-10,10), dns = False):
        #First change figure size
        plt.figure(figsize=(8,6))
        plt.hist(arr, histtype = 'step',linewidth = 1.5, bins=nbins, alpha = 1., range = rng, density = dns)
        plt.xlabel(name)
        plt.ylabel('counts')
        plt.title(name)
        plt.show()

def plot_log(name, arr, nbins=50, rng=(-10,10)):
        #First change figure size
        plt.figure(figsize=(8,6))
        plt.hist(arr, bins=nbins, histtype = 'step', color = 'blue', range=rng)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(name)
        plt.ylabel('counts')
        plt.title(name)
        plt.show()


def plot_cut(name, arr, cut, nbins=50, rng=(-10,10)):
    plt.figure(figsize=(8,6))
    plt.hist(arr, bins=nbins, alpha = 0.5, range = rng)
    plt.hist(arr[cut], bins=nbins, alpha = 0.5, range = rng)
    plt.title(name)
    plt.show()

def double_plot(name1, name2, arr1, arr2,nbins=50, rng=(-10,10)):
    plt.figure(figsize=(20,8))
    plt.subplot(1, 2, 1)
    plt.hist(arr1, bins=nbins, alpha = 1., range = rng)
    plt.ylabel('counts')
    plt.title(name1)
    plt.subplot(1, 2, 2)
    plt.hist(arr2, bins=nbins, alpha = 1., range = rng)
    plt.ylabel('counts')
    plt.title(name2)
    plt.show()

def double_plot_log(name1, name2, arr1, arr2,nbins=50, rng=(-10,10)):
    plt.figure(figsize=(20,8))
    plt.subplot(1, 2, 1)
    plt.hist(arr1, histtype = 'step',linewidth = 1.5, bins=nbins, alpha = 1., range = rng)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('counts')
    plt.title(name1)
    plt.subplot(1, 2, 2)
    plt.hist(arr2, histtype = 'step', linewidth = 1.5, bins=nbins, alpha = 1., range = rng)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('counts')
    plt.title(name2)
    plt.show()

def double_plot_cut(name, arr1, arr2, cut, nbins=50, rng=(-10,10)):
    plt.figure(figsize=(20,8))
    plt.subplot(1, 2, 1)
    plt.hist(arr1, bins=nbins, alpha = 0.5, range = rng)
    plt.hist(arr1[cut],bins=nbins, alpha = 0.5, range = rng)
    plt.ylabel('counts')
    plt.title("muons " + name)
    plt.subplot(1, 2, 2)
    plt.hist(arr2, bins=nbins, alpha = 0.5, range = rng)
    plt.hist(arr2[cut],bins=nbins, alpha = 0.5, range = rng)
    plt.ylabel('counts')
    plt.title("antimuons "+ name)
    plt.show()

def scatter(name, arr):
    plt.figure(figsize=(12,8))
    plt.scatter(arr[:,0],arr[:,1])
    plt.ylabel('y')
    plt.xlabel('x')
    plt.title(name)
    plt.show()
