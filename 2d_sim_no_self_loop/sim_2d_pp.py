# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 08:30:05 2021

@author: t3nemati19
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

dsf=".csv"

N_xy=np.loadtxt("N_xy"+dsf, delimiter=',')
N_x=int(N_xy[0])
N_y=int(N_xy[1])
N_nodes = N_x * N_y

site_indices_inv_array = np.loadtxt('site_indices_inv.csv', delimiter=',', dtype=int)
site_indices = dict({})
site_indices_inv = []
# site_indices[(0,0)]=i
for i in range(np.shape(site_indices_inv_array)[0]):
    n_x = int(site_indices_inv_array[i,0])
    n_y = int(site_indices_inv_array[i,1])
    site_indices[(n_x, n_y)] = i
    site_indices_inv.append((n_x, n_y))
    

iter_data=np.loadtxt("iter_data"+dsf, delimiter=',', dtype=int)
N_blocks=iter_data[0]
N_iter_perblock_pernode=iter_data[1]

fix_iter=np.loadtxt("fix_iter"+dsf, delimiter=',', dtype=int)
fix_time=np.loadtxt("fix_time"+dsf, delimiter=',', dtype=int)
abs_time=np.loadtxt("abs_time"+dsf, delimiter=',', dtype=int)
ext_iter=1-fix_iter
ext_time=ext_iter*abs_time

fix_probs=np.zeros(N_blocks)
fix_times=np.zeros(N_blocks)
abs_times=np.zeros(N_blocks)
ext_times=np.zeros(N_blocks)

local_FP_blocks = np.zeros((N_nodes, N_blocks))
local_FT_blocks = np.zeros((N_nodes, N_blocks))

for block_c in range(N_blocks):
    
    begin_ind=block_c*N_iter_perblock_pernode
    end_ind  =block_c*N_iter_perblock_pernode+N_iter_perblock_pernode
    
    block_fix_iter=fix_iter[:,begin_ind:end_ind]
    block_fix_time=fix_time[:,begin_ind:end_ind]
    block_abs_time=abs_time[:,begin_ind:end_ind]
    block_ext_time=ext_time[:,begin_ind:end_ind]
    
    # fix_probs[block_c]=np.mean(np.mean(block_fix_iter,1))
    # abs_times[block_c]=np.mean(np.mean(block_abs_time,1))
    # fix_times[block_c]=np.mean(block_fix_time[np.nonzero(block_fix_time)])
    # ext_times[block_c]=np.mean(block_ext_time[np.nonzero(block_ext_time)])
    
    fix_probs[block_c]=np.sum(block_fix_iter) / (N_nodes * N_iter_perblock_pernode)
    abs_times[block_c]=np.sum(block_abs_time) / (N_nodes * N_iter_perblock_pernode)
    fix_times[block_c]=np.sum(block_fix_time) / np.sum(1 * (block_fix_time > 0) )
    ext_times[block_c]=np.sum(block_ext_time) / np.sum(1 * (block_ext_time > 0) )
    
    for node_index in range(N_nodes):
        local_FP_blocks[node_index, block_c] = np.sum(block_fix_iter[node_index,:])/len(block_fix_iter[node_index,:])
        local_FT_blocks[node_index, block_c] = np.sum(block_fix_time[node_index,:])/np.sum(1 * (block_fix_iter[node_index,:] > 0) )

np.savetxt("fix_probs_blocks"+dsf, fix_probs, delimiter=',')
np.savetxt("abs_times_blocks"+dsf, abs_times, delimiter=',')
np.savetxt("fix_times_blocks"+dsf, fix_times, delimiter=',')
np.savetxt("ext_times_blocks"+dsf, ext_times, delimiter=',')


# plot local fix prob
local_FP_avg = np.mean(local_FP_blocks, axis=1)
local_FP_err = np.std(local_FP_blocks, axis=1)
# local_FP_err = np.std(local_FP_blocks, axis=1)/np.sqrt(N_blocks)
local_FT_avg = np.mean(local_FT_blocks, axis=1)
local_FT_err = np.std(local_FT_blocks, axis=1)
# local_FT_err = np.std(local_FT_blocks, axis=1)/np.sqrt(N_blocks)

FP_avg_mat_plot = np.zeros((N_x,N_y))
FT_avg_mat_plot = np.zeros((N_x,N_y))

for n_x in range(N_x):
    for n_y in range(N_y):
        index = site_indices[(n_x, n_y)]
        FP_avg_mat_plot[n_x,n_y] = local_FP_avg[index]
        FT_avg_mat_plot[n_x,n_y] = local_FT_avg[index]

FP_norm_avg_std = np.zeros(2)
FP_norm_avg_std[0] = np.mean(fix_probs*N_nodes)
FP_norm_avg_std[1] = np.std(fix_probs*N_nodes)
np.savetxt('FP_norm_avg_std.csv', X=FP_norm_avg_std, delimiter=',')

# plt.figure()
# x = np.linspace(0,N_nodes-1,N_nodes)
# plt.errorbar(x, y=local_FP_avg*N_nodes, yerr=local_FP_err*N_nodes, fmt='o')
# plt.title('FP')
# plt.savefig('FP_vs_index.PNG', dpi=300)
# plt.close()

# plt.figure()
# x = np.linspace(0,N_nodes-1,N_nodes)
# plt.errorbar(x, y=local_FT_avg, yerr=local_FT_err, fmt='o')
# plt.yscale('log')
# plt.title('FT')
# plt.savefig('FT_vs_index.PNG', dpi=300)
# plt.close()



plt.figure()
# plt.imshow(mat, interpolation='nearest', cmap=cm.Greys_r)
plt.imshow(FP_avg_mat_plot*N_nodes, interpolation='nearest', cmap=cm.Blues)
plt.title('FP*N')
plt.colorbar()
plt.axis('equal')
# plt.savefig('config.png', dpi=200)
# plt.show()
plt.savefig('FP_local.PNG', dpi=300)
plt.close()

# plot local fix prob

plt.figure()
# plt.imshow(mat, interpolation='nearest', cmap=cm.Greys_r)
plt.imshow(FT_avg_mat_plot, interpolation='nearest', cmap=cm.Blues)
plt.title('FT')
plt.colorbar()
plt.axis('equal')
# plt.savefig('config.png', dpi=200)
# plt.show()
plt.savefig('FT_local.PNG', dpi=300)
plt.close()
# plot local fix prob
        



