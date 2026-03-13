# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 11:20:36 2023

@author: Nemat002
"""

import numpy as np
from scipy import sparse
from scipy import linalg
import scipy.sparse.linalg as sp_linalg
from time import time

# N_nodes = 9
N_nodes = int(np.loadtxt('N_nodes.txt', delimiter=','))

# rA_vec = np.array([0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.6])
# rA_vec = np.array([0.9, 1.0, 1.1])
rA_vec = np.loadtxt('rA_vec.csv', delimiter=',')
# rA_vec = np.array([0.8, 1.0, 1.2])
np.savetxt('rA_vec.csv', rA_vec, delimiter=',')



line_c =0
fitness_data = np.zeros([1, N_nodes])
mA_vec_tot = []
rA_vec_tot = []

for rA_c in range(len(rA_vec)):
    
    rA = rA_vec[rA_c]
    margin = 0.1
    dA_max = rA - margin
    mA_vec = []
    delta_dA = 0.05
    
    dA = 0
    while(dA < dA_max+(1e-7)):
        mA_vec.append(2*dA/(N_nodes-1.))
        mA_vec_tot.append(2*dA/(N_nodes-1.))
        rA_vec_tot.append(rA)
        dA=dA+delta_dA
    
    for mA_c in range(len(mA_vec)):
        mA = mA_vec[mA_c]
        
        fitness_config = np.zeros([1, N_nodes])
        for node_c in range(N_nodes):
            fitness_config[0,node_c] = rA + mA * (node_c - (N_nodes-1)/2)
        
        if line_c == 0:            
            fitness_data = fitness_config.copy()
        else:
            fitness_data=np.concatenate((fitness_data, fitness_config), axis=0)
        
        line_c+=1

np.savetxt('fitness_data.csv', fitness_data, delimiter=',', fmt='%2.8f\t')
mA_vec_tot = np.array(mA_vec_tot)
np.savetxt('mA_vec_tot.csv', mA_vec_tot, delimiter=',', fmt='%2.8f\t')
rA_vec_tot = np.array(rA_vec_tot)
np.savetxt('rA_vec_tot.csv', rA_vec_tot, delimiter=',', fmt='%2.8f\t')