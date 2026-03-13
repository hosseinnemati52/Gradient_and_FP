# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 09:32:05 2021

@author: t3nemati19
"""

import numpy as np



def neighborhood_matrix(config_keyword, N):
    if config_keyword=="cycle":
        matrix=np.zeros((N,N))
        for i in range(N):
            matrix[i,(i-1)%N]=1
            matrix[i,(i+1)%N]=1
    
    for i in range(N):
        matrix[i,:]/=np.sum(matrix[i,:])
    return matrix

def neighborhood_matrix_2d(config_keyword, N_xy, BC_key):
    N_x=N_xy[0]
    N_y=N_xy[1]
    N_tot=N_x*N_y
    if config_keyword=="sq_lattice":
        
        if BC_key=='open':
            # matrix=np.zeros((N_tot,N_tot))
            
            # site_indices=dict({})
            # i=0
            # site_indices[(0,0)]=i
            # for n_y in range(N_y):
            #     for n_x in range(N_x):
            #         site_indices[(n_x, n_y)]=i
                    
            #         if n_y==0:
            #             if n_x==0:
            #                 i_r=i+1
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_d]=1
                            
            #             elif n_x==N_x-1:
            #                 i_l=i-1
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_d]=1
                            
            #             else:
            #                 i_r=i+1
            #                 i_l=i-1
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_d]=1
                            
            #         elif n_y==N_y-1:
            #             if n_x==0:
                            
            #                 i_r=i+1
            #                 i_u=i-N_x
                            
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_u]=1
                            
            #             elif n_x==N_x-1:
                            
            #                 i_l=i-1
            #                 i_u=i-N_x
                            
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_u]=1
                            
            #             else:
                            
            #                 i_l=i-1
            #                 i_r=i+1
            #                 i_u=i-N_x
                            
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_u]=1
                            
            #         else:
            #             if n_x==0:
                            
            #                 i_r=i+1
            #                 i_u=i-N_x
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_u]=1
            #                 matrix[i,i_d]=1
                            
            #             elif n_x==N_x-1:
                            
            #                 i_l=i-1
            #                 i_u=i-N_x
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_u]=1
            #                 matrix[i,i_d]=1
                            
            #             else:
                            
            #                 i_l=i-1
            #                 i_r=i+1
            #                 i_u=i-N_x
            #                 i_d=i+N_x
                            
            #                 matrix[i,i_l]=1
            #                 matrix[i,i_r]=1
            #                 matrix[i,i_u]=1
            #                 matrix[i,i_d]=1
                    
            #         i+=1
            pass
                
        elif BC_key=='periodic':
            
            matrix=np.zeros((N_tot,N_tot))
            
            site_indices     = dict({})
            site_indices_inv = []
            i=0
            # site_indices[(0,0)]=i
            for n_x in range(N_x):
                for n_y in range(N_y):
                    site_indices[(n_x, n_y)]=i
                    site_indices_inv.append((n_x, n_y))
                    i+=1
            
            for n_x in range(N_x):
                for n_y in range(N_y):
                    
                    i=site_indices[(n_x, n_y)]
                    i_d=site_indices[((n_x+1)%N_x, n_y)]
                    i_u=site_indices[((n_x-1)%N_x, n_y)]
                    i_l=site_indices[(n_x, (n_y-1)%N_y)]
                    i_r=site_indices[(n_x, (n_y+1)%N_y)]
                    
                    matrix[i,i_r]=1
                    matrix[i,i_l]=1
                    matrix[i,i_u]=1
                    matrix[i,i_d]=1
                    
            
            for self_link_c in range(N_tot):
                matrix[self_link_c,self_link_c]=0
        
    elif config_keyword=="cylinder":
        
        if BC_key=='x_per__y_open':
            
            matrix=np.zeros((N_tot,N_tot))
            
            site_indices     = dict({})
            site_indices_inv = []
            i=0
            # site_indices[(0,0)]=i
            for n_x in range(N_x):
                for n_y in range(N_y):
                    site_indices[(n_x, n_y)]=i
                    site_indices_inv.append((n_x, n_y))
                    i+=1
            
            for n_x in range(N_x):
                for n_y in range(N_y):
                    
                    i=site_indices[(n_x, n_y)]
                    i_d=site_indices[((n_x+1)%N_x, n_y)]
                    i_u=site_indices[((n_x-1)%N_x, n_y)]
                    matrix[i,i_u]=1
                    matrix[i,i_d]=1
                    
                    if n_y>0:
                        # i_l=site_indices[(n_x, (n_y-1)%N_y)]
                        i_l=site_indices[(n_x, n_y-1)]
                        matrix[i,i_l]=1
                    else:
                        pass
                    
                    if n_y<N_y-1:
                        # i_r=site_indices[(n_x, (n_y+1)%N_y)]
                        i_r=site_indices[(n_x, n_y+1)]
                        matrix[i,i_r]=1
                    else:
                        pass
            
            for self_link_c in range(N_tot):
                matrix[self_link_c,self_link_c]=0
                    
    
    for i in range(N_tot):
        matrix[i,:]=matrix[i,:]/np.sum(matrix[i,:])
    
    
    return matrix , site_indices, site_indices_inv

def neighborhood_lists_weights(matrix):
    
    lists=  dict({})
    weights=dict({})
    N=np.shape(matrix)[0]
    for i in range(N):
        lists[i]=[]
        weights[i]=[]
        for j in range(N):
            if matrix[i,j]>0:
                lists[i].append(j)
                weights[i].append(matrix[i,j])
    return (lists,weights)

def write_geom(config_keyword,N):
    config_file = open("config"+dsf, "wt")
    config_file.write(config_keyword)
    config_file.close()
    np.savetxt("N"+dsf, np.array([N]), fmt='%d')
    return

def write_geom_2d(config_keyword, N_xy, BC_key):
    config_file = open("config"+dsf, "wt")
    config_file.write(config_keyword)
    config_file.close()
    
    BC_file = open("BC"+dsf, "wt")
    BC_file.write(BC_key)
    BC_file.close()
    
    # np.savetxt("N_xy"+dsf, np.array(N_xy), fmt='%d')
    return

def fitness_setup(N):
    fitnesses={'resident':np.zeros(N), 'mutant':np.zeros(N)}
    
    r_res=1
    r_mut=1.1
    np.savetxt("r_res_r_mut"+dsf, np.array([r_res,r_mut]), fmt='%.3e', delimiter=',')
    
    s_res=0.5
    s_mut=0.5
    np.savetxt("s_res_s_mut"+dsf, np.array([s_res,s_mut]), fmt='%.3e', delimiter=',')
    
    #Periodic fitnesses
    T=8
    np.savetxt("period"+dsf, np.array([T]), fmt='%d', delimiter=',')
    N_period=int(N/T)
    i=0
    for p_ind in range(N_period):
        for i_dum in range(int(T/2)):
            fitnesses['resident'][i]=r_res+s_res
            fitnesses['mutant'][i]=r_mut+s_mut
            i=i+1
        for i_dum in range(int(T/2)):
            fitnesses['resident'][i]=r_res-s_res
            fitnesses['mutant'][i]=r_mut-s_mut
            i=i+1
        
    fitness_table=np.zeros((2,N))
    fitness_table[0,:]=fitnesses['resident']
    fitness_table[1,:]=fitnesses['mutant']
    return fitness_table
    
def fitness_setup_2d(N_xy, site_indices, site_indices_inv):
    N_x=N_xy[0]
    N_y=N_xy[1]
    
    N_tot=N_x*N_y
    
    fitnesses={'resident':np.zeros(N_tot), 'mutant':np.zeros(N_tot)}
    
    r_res=1
    r_mut=1.5
    np.savetxt("r_res_r_mut"+dsf, np.array([r_res,r_mut]), fmt='%.3e', delimiter=',')
    
    sigma = 0.6
    
    scenario = np.loadtxt('scenario.csv', dtype=int)
    if scenario==1:
        s_res = sigma
        s_mut = sigma
    elif scenario==2:
        s_res = 0
        s_mut = sigma
    elif scenario==3:
        s_res = sigma
        s_mut = 0
    elif scenario==4:
        s_res =  sigma
        s_mut = -sigma
        
    np.savetxt("s_res_s_mut"+dsf, np.array([s_res,s_mut]), fmt='%.3e', delimiter=',')
    
    #Green and Red nodes (Green: r+s, Red: r-s)
    # site_colors=dict({})
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if n_x<=1:
    #             site_colors[(n_x, n_y)]='g'
    #         else:
    #             site_colors[(n_x, n_y)]='r'
    
    
    # #Periodic colors
    # T=2
    # np.savetxt("period"+dsf, np.array([T]), fmt='%d', delimiter=',')
    # N_period=int(N_x/T)
    # site_colors=dict({})
    # i=0
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if np.floor(n_x/int(T/2))%2==0:
    #             site_colors[(n_x, n_y)]='g'
    #             i=i+1
    #         if np.floor(n_x/int(T/2))%2==1:
    #             site_colors[(n_x, n_y)]='r'
    #             i=i+1
    
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if site_colors[(n_x, n_y)]=='g':
    #             fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res+s_res
    #             fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut+s_mut
    #         elif site_colors[(n_x, n_y)]=='r':
    #             fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res-s_res
    #             fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut-s_mut
    
    # random config
    site_colors=dict({})
    try:
        shuffled_indices = np.loadtxt('config_g_r.csv', dtype = int, delimiter=',')
    except:
        shuffled_indices = np.array(list(range(N_tot)))
        np.random.shuffle(shuffled_indices)
        np.savetxt('config_g_r.csv', shuffled_indices, fmt="%d", delimiter=',')
        
    i=0
    
    for node_c in range(int(N_tot/2)):
        site_colors[site_indices_inv[shuffled_indices[node_c]]] = 'g'
    for node_c in range(int(N_tot/2), N_tot):
        site_colors[site_indices_inv[shuffled_indices[node_c]]] = 'r'
        
    
    for n_x in range(N_x):
        for n_y in range(N_y):
            if site_colors[(n_x, n_y)]=='g':
                fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res+s_res
                fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut+s_mut
            elif site_colors[(n_x, n_y)]=='r':
                fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res-s_res
                fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut-s_mut
        
    fitness_table=np.zeros((2,N_tot))
    fitness_table[0,:]=fitnesses['resident']
    fitness_table[1,:]=fitnesses['mutant']
    return fitness_table

def fitness_setup_cyl_grad(N_xy, site_indices, site_indices_inv):
    N_x=N_xy[0]
    N_y=N_xy[1]
    
    N_tot=N_x*N_y
    
    fitnesses={'resident':np.zeros(N_tot), 'mutant':np.zeros(N_tot)}
    
    # r_res=1
    # r_mut=1.5
    # np.savetxt("r_res_r_mut"+dsf, np.array([r_res,r_mut]), fmt='%.3e', delimiter=',')
    
    # sigma = 0.6
    
    r_s_data = np.loadtxt('r_s_data.csv', dtype=float, delimiter=',')
    r_res = r_s_data[0]
    r_mut = r_s_data[1]
    sigma = r_s_data[2]
    
    scenario = np.loadtxt('scenario.csv', dtype=int)
    if scenario==1:
        for n_x in range(N_x):
            for n_y in range(N_y):
                fitnesses['resident'][site_indices[(n_x, n_y)]] = (r_res-sigma) + 2*sigma*n_y/(N_y-1)
                fitnesses['mutant'][site_indices[(n_x, n_y)]]   = (r_mut-sigma) + 2*sigma*n_y/(N_y-1)
    elif scenario==2:
        for n_x in range(N_x):
            for n_y in range(N_y):
                fitnesses['resident'][site_indices[(n_x, n_y)]] =  r_res
                fitnesses['mutant'][site_indices[(n_x, n_y)]]   = (r_mut-sigma) + 2*sigma*n_y/(N_y-1)
    elif scenario==3:
        for n_x in range(N_x):
            for n_y in range(N_y):
                fitnesses['resident'][site_indices[(n_x, n_y)]] = (r_res-sigma) + 2*sigma*n_y/(N_y-1)
                fitnesses['mutant'][site_indices[(n_x, n_y)]]   =  r_mut
    elif scenario==4:
        for n_x in range(N_x):
            for n_y in range(N_y):
                fitnesses['resident'][site_indices[(n_x, n_y)]] = (r_res-sigma) + 2*sigma*n_y/(N_y-1)
                fitnesses['mutant'][site_indices[(n_x, n_y)]]   = (r_mut+sigma) - 2*sigma*n_y/(N_y-1)
        
    # np.savetxt("s_res_s_mut"+dsf, np.array([s_res,s_mut]), fmt='%.3e', delimiter=',')
    
    #Green and Red nodes (Green: r+s, Red: r-s)
    # site_colors=dict({})
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if n_x<=1:
    #             site_colors[(n_x, n_y)]='g'
    #         else:
    #             site_colors[(n_x, n_y)]='r'
    
    
    # #Periodic colors
    # T=2
    # np.savetxt("period"+dsf, np.array([T]), fmt='%d', delimiter=',')
    # N_period=int(N_x/T)
    # site_colors=dict({})
    # i=0
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if np.floor(n_x/int(T/2))%2==0:
    #             site_colors[(n_x, n_y)]='g'
    #             i=i+1
    #         if np.floor(n_x/int(T/2))%2==1:
    #             site_colors[(n_x, n_y)]='r'
    #             i=i+1
    
    # for n_y in range(N_y):
    #     for n_x in range(N_x):
    #         if site_colors[(n_x, n_y)]=='g':
    #             fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res+s_res
    #             fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut+s_mut
    #         elif site_colors[(n_x, n_y)]=='r':
    #             fitnesses['resident'][site_indices[(n_x, n_y)]]=r_res-s_res
    #             fitnesses['mutant'][site_indices[(n_x, n_y)]]=r_mut-s_mut
    
    # random config
    
    
        
    fitness_table=np.zeros((2,N_tot))
    fitness_table[0,:]=fitnesses['resident']
    fitness_table[1,:]=fitnesses['mutant']
    return fitness_table


def write_iteration_data(N_blocks,N_iter_perblock_pernode):
    np.savetxt("iter_data"+dsf, np.array([N_blocks,N_iter_perblock_pernode]), fmt='%d', delimiter=',')
    return

def solver(N,neighbors_mat,neighbors_lists,neighbors_weights,fitness_table,N_blocks,N_iter_perblock_pernode):
    
    try:
        saving_key=np.loadtxt("saving_key"+dsf,dtype=int, delimiter=',')
        
        fix_iter=np.loadtxt("fix_iter"+dsf,dtype=int, delimiter=',')
        fix_time=np.loadtxt("fix_time"+dsf,dtype=int, delimiter=',')
        abs_time=np.loadtxt("abs_time"+dsf,dtype=int, delimiter=',')
    except OSError:
        saving_key=np.zeros((N,N_blocks), dtype=int)
        np.savetxt("saving_key"+dsf,saving_key,fmt="%d", delimiter=',')
        
        fix_iter=np.zeros((N,N_blocks*N_iter_perblock_pernode))
        fix_time=np.zeros((N,N_blocks*N_iter_perblock_pernode))
        abs_time=np.zeros((N,N_blocks*N_iter_perblock_pernode))
        np.savetxt("fix_iter"+dsf,fix_iter,fmt="%d", delimiter=',')
        np.savetxt("fix_time"+dsf,fix_time,fmt="%d", delimiter=',')
        np.savetxt("abs_time"+dsf,abs_time,fmt="%d", delimiter=',')
        
    # saving_key=np.zeros((N,N_blocks), dtype=int)
    # np.savetxt("saving_key.txt",saving_key,fmt="%d")
    
    # fix_iter=np.zeros((N,N_blocks*N_iter_perblock_pernode))
    # fix_time=np.zeros((N,N_blocks*N_iter_perblock_pernode))
    # abs_time=np.zeros((N,N_blocks*N_iter_perblock_pernode))
    # np.savetxt("fix_iter.txt",fix_iter,fmt="%d")
    # np.savetxt("fix_time.txt",fix_time,fmt="%d")
    # np.savetxt("abs_time.txt",abs_time,fmt="%d")
    
    
    for block_c in range(N_blocks):
    
        for node_c in range(N):
            
            print('block_c: '+str(block_c+1)+'/'+str(N_blocks)+' , node_c: '+str(node_c)+'/'+str(N-1))
            
            if saving_key[node_c,block_c]==0:
                
                for iter_c in range(N_iter_perblock_pernode):
                    
                    # print('block_c: '+str(block_c+1)+'/'+str(N_blocks)+' , node_c: '+str(node_c)+'/'+str(N-1)+', iter: '+str(iter_c+1)+'/'+str(N_iter_perblock_pernode))
                    
                    #initialization
                    ResMut=np.zeros(N,dtype=int)
                    fitness=(fitness_table[0,:]).copy()
                    
                    ResMut[node_c]=1
                    fitness[node_c]=(fitness_table[ResMut[node_c],node_c]).copy()
                    t=0
                    
                    while np.std(ResMut)>0:
                        
                        probs=(fitness.copy())/np.sum(fitness.copy())
                        birth_node=np.random.choice(N,p=probs)
                        # death_node=np.random.choice(neighbors_lists[birth_node],p=neighbors_weights[birth_node])
                        # death_node=neighbors_lists[birth_node][np.random.randint(4)]
                        death_node=neighbors_lists[birth_node][np.random.randint(len(neighbors_lists[birth_node]))]
                        
                        ResMut[death_node]=ResMut[birth_node].copy()
                        fitness[death_node]=(fitness_table[ResMut[death_node],death_node]).copy()
                        t+=1
                        
                    row_ind=node_c
                    col_ind=block_c*N_iter_perblock_pernode+iter_c
                    
                    abs_time[row_ind,col_ind]=t
                    if np.mean(ResMut)>0.5:
                        fix_iter[row_ind,col_ind]=1
                        fix_time[row_ind,col_ind]=t
                
                np.savetxt("fix_iter"+dsf,fix_iter,fmt="%d", delimiter=',')
                np.savetxt("fix_time"+dsf,fix_time,fmt="%d", delimiter=',')
                np.savetxt("abs_time"+dsf,abs_time,fmt="%d", delimiter=',')
                
                saving_key[node_c,block_c]=1
                np.savetxt("saving_key"+dsf,saving_key,fmt="%d", delimiter=',')
    
    return (fix_iter,fix_time,abs_time)

dsf=".csv"
# config_keyword="cycle"
# config_keyword="sq_lattice"
config_keyword="cylinder"


if config_keyword=="cycle":
    N=16
    write_geom(config_keyword,N)
    neighbors_mat=neighborhood_matrix(config_keyword,N)
elif config_keyword=="sq_lattice":
    N_x=8
    N_y=8
    N_tot=N_x*N_y
    N_xy=[N_x, N_y]
    # BC_key='open'
    BC_key='periodic'
elif config_keyword=="cylinder":
    
    N_xy_read = np.loadtxt('N_xy.csv', delimiter=',', dtype=int)
    
    N_x = int(N_xy_read[0])
    N_y = int(N_xy_read[1])
    N_tot=N_x*N_y
    N_xy=[N_x, N_y]
    # BC_key='open'
    BC_key='x_per__y_open'
    
    write_geom_2d(config_keyword, N_xy, BC_key)
    
    neighbors_mat, site_indices, site_indices_inv = neighborhood_matrix_2d(config_keyword, N_xy, BC_key)
    np.savetxt('site_indices_inv.csv', X=np.array(site_indices_inv), fmt='%d', delimiter=',')
    
    # if BC_key=='open':
    #     pass
    # elif BC_key=='periodic':
    #     neighbors_mat=periodic_BC(N_xy,neighbors_mat,site_indices)
    
    

(neighbors_lists,neighbors_weights)=neighborhood_lists_weights(neighbors_mat)



if config_keyword=="cycle":
    fitness_table=fitness_setup(N)
    
elif config_keyword=="sq_lattice":
    fitness_table=fitness_setup_2d(N_xy, site_indices, site_indices_inv)

elif config_keyword=="cylinder":
    fitness_table=fitness_setup_cyl_grad(N_xy, site_indices, site_indices_inv)
    


N_blocks=10
N_iter_perblock_pernode=400
write_iteration_data(N_blocks,N_iter_perblock_pernode)


if config_keyword=="cycle":
    (fix_iter,fix_time,abs_time)=solver(N,neighbors_mat,neighbors_lists,neighbors_weights,fitness_table,N_blocks,N_iter_perblock_pernode)
elif config_keyword=="sq_lattice":
    (fix_iter,fix_time,abs_time)=solver(N_tot,neighbors_mat,neighbors_lists,neighbors_weights,fitness_table,N_blocks,N_iter_perblock_pernode)
elif config_keyword=="cylinder":
    (fix_iter,fix_time,abs_time)=solver(N_tot,neighbors_mat,neighbors_lists,neighbors_weights,fitness_table,N_blocks,N_iter_perblock_pernode)

# # import sim_pp
# exec(open('sim_pp.py').read()) 

# import sim_2d_pp
exec(open('sim_2d_pp.py').read()) 




            