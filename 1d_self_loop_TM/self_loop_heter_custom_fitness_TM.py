#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 19:37:17 2021

@author: hossein
"""

# from numba import jit
import numpy as np
from scipy import sparse
from scipy import linalg
import scipy.sparse.linalg as sp_linalg
from time import time
import seaborn as sns
import matplotlib.pyplot as plt
import os



# def Caley_Hamilton_inv(A):
#     d=3
#     return A_inv
    

# @jit
def Q_R_maker(fitnesses_data):


    # fitnesses_data=np.array([ [1,1,1] , [2,3,4] ]) #row[0]: resident, row[1]: mutant
    
    
    N_nodes = len(fitnesses_data[0,:])
    
    L_Q = N_nodes*(N_nodes-1)
    
    #R Calculation
    R=np.zeros((L_Q,2), np.longdouble) #there are two abrorbing states
    
    for i in range(N_nodes):
 
        F_tot=np.sum(fitnesses_data[0,:])-fitnesses_data[0,i]+fitnesses_data[1,i]
        i_p=(i+1)%N_nodes
        i_m=(i-1)%N_nodes
        R[i,0]=(fitnesses_data[0,i_m]+fitnesses_data[0,i_p])/(2*F_tot)


        i_p=(i+1)%N_nodes
        i_m=(i-1)%N_nodes        
        i_2m=(i-2)%N_nodes
        F_tot=np.sum(fitnesses_data[1,:])-fitnesses_data[1,i_m]+fitnesses_data[0,i_m]
        R[N_nodes*(N_nodes-1)-N_nodes+i,1]=(fitnesses_data[1,i]+fitnesses_data[1,i_2m])/(2*F_tot)
    
    
    
    
    
    #Q Calculation
    Q=np.zeros((L_Q,L_Q), np.longdouble)
    for block_c_y in range(N_nodes-1):
        
        ###########################################################################
        # diagonal blocks
        subQ=np.zeros((N_nodes,N_nodes))
        num_A=block_c_y+1 # num_A is number of mutant elements
        num_B=N_nodes-num_A
        # A_chain_begin_pos is beginnging position of the mutant chain (0<=A_chain_begin_pos<=N_nodes-1)
        for A_chain_begin_pos in range(N_nodes):
    
            cycle=np.zeros(N_nodes,dtype='int')
            fitness=np.zeros(N_nodes)
            for A_pos_index in range(num_A):
                j=(A_chain_begin_pos+A_pos_index)%N_nodes
                cycle[j]=1
                fitness[j]=fitnesses_data[1,j]
            A_chain_end_pos=(A_chain_begin_pos+num_A-1)%N_nodes
            B_chain_begin_pos=(A_chain_end_pos+1)%N_nodes
            for B_pos_index in range(num_B):
                j=(B_chain_begin_pos+B_pos_index)%N_nodes
                cycle[j]=0
                fitness[j]=fitnesses_data[0,j]
            F_tot=np.sum(fitness)    
            prob_numenator=0
            for j in range(N_nodes):
                j_p=(j+1)%N_nodes
                j_m=(j-1)%N_nodes
                prob_numenator+=fitness[j]*(cycle[j_p]+cycle[j_m])*cycle[j]
                prob_numenator+=fitness[j]*(2-cycle[j_p]-cycle[j_m])*(1-cycle[j])            
            prob=prob_numenator/(2*F_tot)            
            subQ[A_chain_begin_pos,A_chain_begin_pos]=prob
        
        insert_subQ_y_begin_index=block_c_y*N_nodes
        insert_subQ_y_end_index=block_c_y*N_nodes+N_nodes
        insert_subQ_x_begin_index=insert_subQ_y_begin_index  #it is a diagonal block
        insert_subQ_x_end_index=insert_subQ_y_end_index      #it is a diagonal block
        
        Q[insert_subQ_y_begin_index:insert_subQ_y_end_index,insert_subQ_x_begin_index:insert_subQ_x_end_index]=subQ
        #################################################################################
        
        #################################################################################
        #Right to diagonal blocks (add an "A" to the "A" chain)
        if block_c_y<N_nodes-2: # if it is not the last row of blocks do the whole thing
            subQ=np.zeros((N_nodes,N_nodes))
            for A_chain_begin_pos in range(N_nodes):
                k=A_chain_begin_pos
                k_m=(k-1)%N_nodes
                
                cycle=np.zeros(N_nodes,dtype='int')
                fitness=np.zeros(N_nodes)
                for A_pos_index in range(num_A):
                    j=(A_chain_begin_pos+A_pos_index)%N_nodes
                    cycle[j]=1
                    fitness[j]=fitnesses_data[1,j]
                A_chain_end_pos=(A_chain_begin_pos+num_A-1)%N_nodes
                B_chain_begin_pos=(A_chain_end_pos+1)%N_nodes
                for B_pos_index in range(num_B):
                    j=(B_chain_begin_pos+B_pos_index)%N_nodes
                    cycle[j]=0
                    fitness[j]=fitnesses_data[0,j]                    
                prob_numenator_diag=0
                prob_numenator_left_diag=0
                for j in range(N_nodes):
                    j_p=(j+1)%N_nodes
                    j_m=(j-1)%N_nodes
                    prob_numenator_diag+=fitness[j]*(1-cycle[j_p])*cycle[j]
                    prob_numenator_left_diag+=fitness[j]*(1-cycle[j_m])*cycle[j]                    
                F_tot=np.sum(fitness)
                prob_diag=prob_numenator_diag/(2*F_tot) #diagonal elements of subQ
                prob_left_diag=prob_numenator_left_diag/(2*F_tot) #left of diagonal elements of subQ
                
                subQ[k,k]=prob_diag
                subQ[k,k_m]=prob_left_diag
            
            insert_subQ_y_begin_index=block_c_y*N_nodes
            insert_subQ_y_end_index=block_c_y*N_nodes+N_nodes
            insert_subQ_x_begin_index=insert_subQ_y_begin_index+N_nodes  #it is a right of diagonal block
            insert_subQ_x_end_index=insert_subQ_y_end_index+N_nodes      #it is a right of diagonal block
            
            Q[insert_subQ_y_begin_index:insert_subQ_y_end_index,insert_subQ_x_begin_index:insert_subQ_x_end_index]=subQ
        #########################################################################################    
        
        ########################################################################################
        #Left to diagonal blocks (remove an "A" from the "A" chain)
        if block_c_y>0: # if it is not the first row of blocks do the whole thing
            subQ=np.zeros((N_nodes,N_nodes))
            for A_chain_begin_pos in range(N_nodes):                
                k=A_chain_begin_pos
                k_p=(k+1)%N_nodes
                cycle=np.zeros(N_nodes,dtype='int')
                fitness=np.zeros(N_nodes)                
                for A_pos_index in range(num_A):
                    j=(A_chain_begin_pos+A_pos_index)%N_nodes
                    cycle[j]=1
                    fitness[j]=fitnesses_data[1,j]                
                A_chain_end_pos=(A_chain_begin_pos+num_A-1)%N_nodes
                B_chain_begin_pos=(A_chain_end_pos+1)%N_nodes
                for B_pos_index in range(num_B):
                    j=(B_chain_begin_pos+B_pos_index)%N_nodes
                    cycle[j]=0
                    fitness[j]=fitnesses_data[0,j]                
                prob_numenator_diag=0
                prob_numenator_right_diag=0    
                for j in range(N_nodes):
                    j_p=(j+1)%N_nodes
                    j_m=(j-1)%N_nodes
                    prob_numenator_diag+=fitness[j]*cycle[j_m]*(1-cycle[j])
                    prob_numenator_right_diag+=fitness[j]*cycle[j_p]*(1-cycle[j])
                F_tot=np.sum(fitness)
                prob_diag=prob_numenator_diag/(2*F_tot) #diagonal elements of subQ
                prob_right_diag=prob_numenator_right_diag/(2*F_tot) #right of diagonal elements of subQ
                
                subQ[k,k]=prob_diag
                subQ[k,k_p]=prob_right_diag
            
            insert_subQ_y_begin_index=block_c_y*N_nodes
            insert_subQ_y_end_index=block_c_y*N_nodes+N_nodes
            insert_subQ_x_begin_index=insert_subQ_y_begin_index-N_nodes  #it is a left of diagonal block
            insert_subQ_x_end_index=insert_subQ_y_end_index-N_nodes      #it is a left of diagonal block
            
            Q[insert_subQ_y_begin_index:insert_subQ_y_end_index,insert_subQ_x_begin_index:insert_subQ_x_end_index]=subQ
        #########################################################################################    
        
    return (Q,R)

def indexFinder(N_nodes, length, start):
    index = -1
    for m_length in range(1, N_nodes):
        for s_point in range(0, N_nodes - m_length + 1):
            index += 1
            if ((m_length == length) and (s_point == start)):
                
                return index
            
    return np.nan
    

def Q_R_maker_self_loop(fitnesses_data):


    # fitnesses_data=np.array([ [1,1,1] , [2,3,4] ]) #row[0]: resident, row[1]: mutant
    
    
    N_nodes = len(fitnesses_data[0,:])
    
    # L_Q = N_nodes*(N_nodes-1)
    L_Q = int(0.5*N_nodes*(N_nodes+1)-1)
    
    #R Calculation
    R=np.zeros((L_Q,2), np.longdouble) #there are two abrorbing states
    
    # extinction    
    for i in range(N_nodes):

        F_tot=np.sum(fitnesses_data[0,:])-fitnesses_data[0,i]+fitnesses_data[1,i]
        if i==0:
            i_p=(i+1)
            R[i,0]=(fitnesses_data[0,i_p])/(2*F_tot)
        elif i==(N_nodes-1):
            i_m=(i-1)
            R[i,0]=(fitnesses_data[0,i_m])/(2*F_tot)
        else:
            i_p=(i+1)
            i_m=(i-1)
            R[i,0]=(fitnesses_data[0,i_m]+fitnesses_data[0,i_p])/(2*F_tot)
    # extinction
        
    # fixation
    F_tot=np.sum(fitnesses_data[1,:])-fitnesses_data[1,-1]+fitnesses_data[0,-1]
    R[L_Q-2,1]=(fitnesses_data[1,-2])/(2*F_tot)
    
    F_tot=np.sum(fitnesses_data[1,:])-fitnesses_data[1,0]+fitnesses_data[0,0]
    R[L_Q-1,1]=(fitnesses_data[1,1])/(2*F_tot)
    # fixation
        
    
    #Q Calculation
    Q=np.zeros((L_Q,L_Q), np.longdouble)
    
    change_switches =dict()
    change_switches['add_r']=1
    change_switches['del_r']=1
    change_switches['add_l']=1
    change_switches['del_l']=1
    
    change_probs =dict()
    change_probs['add_r']=0.
    change_probs['del_r']=0.
    change_probs['add_l']=0.
    change_probs['del_l']=0.
    
    indFrom = -1
    for mutant_length in range(1, N_nodes):
        for start_point in range(0, N_nodes - mutant_length + 1):
            
            indFrom += 1
            end_point = start_point + mutant_length - 1 #(included)
            
            config = np.zeros(N_nodes, dtype=int)
            config[start_point:end_point+1] = 1
            
            sum_non_self = 0
            
            ## Switches
            change_switches['add_r']=1
            change_switches['del_r']=1
            change_switches['add_l']=1
            change_switches['del_l']=1        
            
            if (end_point==N_nodes-1) or (mutant_length==N_nodes-1):
                change_switches['add_r']=0
            else:
                change_switches['add_r']=1
            
            if (end_point==N_nodes-1) or (mutant_length==1):
                change_switches['del_r']=0
            else:
                change_switches['del_r']=1
            
            if (start_point==0) or (mutant_length==N_nodes-1):
                change_switches['add_l']=0
            else:
                change_switches['add_l']=1
                
            if (start_point==0) or (mutant_length==1):
                change_switches['del_l']=0
            else:
                change_switches['del_l']=1
            ## Switches
            
            ## Probs
            change_probs['add_r']=0.
            change_probs['del_r']=0.
            change_probs['add_l']=0.
            change_probs['del_l']=0.
            
            F_tot = 0.
            for i in range(N_nodes):
                F_tot += fitnesses_data[config[i],i]
            
            if (change_switches['add_r']>0):
                indTo = indexFinder(N_nodes, mutant_length+1, start_point)
                change_probs['add_r'] = fitnesses_data[1,end_point]/(2*F_tot)
                Q[indFrom, indTo] = change_probs['add_r']
            if (change_switches['del_r']>0):
                indTo = indexFinder(N_nodes, mutant_length-1, start_point)
                change_probs['del_r'] = fitnesses_data[0,end_point+1]/(2*F_tot)
                Q[indFrom, indTo] = change_probs['del_r']
            if (change_switches['add_l']>0):
                indTo = indexFinder(N_nodes, mutant_length+1, start_point-1)
                change_probs['add_l'] = fitnesses_data[1,start_point]/(2*F_tot)
                Q[indFrom, indTo] = change_probs['add_l']
            if (change_switches['del_l']>0):
                indTo = indexFinder(N_nodes, mutant_length-1, start_point+1)
                change_probs['del_l'] = fitnesses_data[0,start_point-1]/(2*F_tot)
                Q[indFrom, indTo] = change_probs['del_l']
            
            Q[indFrom, indFrom] = 1 - change_probs['add_r'] \
                                    - change_probs['del_r'] \
                                    - change_probs['add_l'] \
                                    - change_probs['del_l'] \
                                    - R[indFrom,0] \
                                    - R[indFrom,1]
            ## Probs
        
    return (Q,R)

    
def Q_R_normalizer(Q,R):
    L_Q=np.shape(Q)[0]
    for i in range(L_Q):
        norm=np.sum(Q[i,:])+np.sum(R[i,:])
        Q[i,:]/=norm
        R[i,:]/=norm
    return(Q,R)
    
# @jit
def prob_and_time_calc(Q,R):
    
    start_time=time()
    
    L_Q = np.shape(Q)[0]
    
    # #Numpy linalg
    # Q=np.float32(Q) # if using np linalg
    # # F=np.linalg.inv(np.identity(L_Q)-Q)
    # F=np.linalg.solve(np.identity(L_Q)-Q,np.identity(L_Q))
    # Phi=np.matmul(F, R)
    # Tau=np.sum(F,1)
    # Cond_Tau=np.matmul(F,Phi)/Phi
    
    
    # #Scipy linalg
    # Q=np.longdouble(Q)
    # F=linalg.inv(np.identity(L_Q)-Q)
    # # F=linalg.solve(np.identity(L_Q)-Q,np.identity(L_Q))
    # Phi=F.dot(R)
    # Tau=np.sum(F,1)
    # Cond_Tau=(F.dot(Phi))/Phi
    
    #Scipy sparse linalg
    Q_sp=sparse.csc_matrix(Q)
    R_sp=sparse.csc_matrix(R)
    F_sp=sp_linalg.inv(sparse.csc_matrix(np.identity(L_Q))-Q_sp)
    # F_sp=Caley_Hamilton_inv(sparse.csc_matrix(np.identity(L_Q))-Q_sp)
    # F=linalg.solve(np.identity(L_Q)-Q,np.identity(L_Q))
    
    gg=F_sp.dot(sparse.csc_matrix(np.identity(L_Q))-Q_sp)
    hh=sparse.csc_matrix(np.identity(L_Q))-gg
    print(np.max(hh))
    print(np.min(hh))
    if (np.max(hh)>(1e-6) or np.min(hh)<(-1e-6)):
        print('############################')
        print('more than 1e-6')
        print('exit')
        print('############################')
        os.exit()
    
    Phi_sp=F_sp.dot(R_sp)
    Phi=Phi_sp.todense()
    F=F_sp.todense()
    Cond_Tau=(F.dot(Phi))/Phi
    Tau=np.sum(F,1)
    
    end_time=time()
    print(end_time-start_time)
    
    
    
    
    return (Phi, Tau, Cond_Tau)
    

def scenario_non_0_func():

    # FIX_PROB_A=np.zeros((len(configs), N_nodes))
    # FIX_PROB_B=np.zeros((len(configs), N_nodes))
    # ABS_TIME=  np.zeros((len(configs), N_nodes))
    # FIX_TIME_A=np.zeros((len(configs), N_nodes))
    # FIX_TIME_B=np.zeros((len(configs), N_nodes))
    # EXT_TIME_A=np.zeros((len(configs), N_nodes))
    # EXT_TIME_B=np.zeros((len(configs), N_nodes))
    
    fix_prob = np.zeros((len(configs), N_nodes))
    abs_time = np.zeros((len(configs), N_nodes))
    fix_time = np.zeros((len(configs), N_nodes))
    
    fix_prob_avg = np.zeros((len(configs), 1))
    abs_time_avg = np.zeros((len(configs), 1))
    fix_time_avg = np.zeros((len(configs), 1))
    
    for config_c in range(len(configs)):
        
        fitness_A = fitness_data[config_c, :]
        fitness_B = np.ones(len(fitness_A))
        
        fitnesses_data=np.vstack([fitness_B,fitness_A]) #row[0]: resident, row[1]: mutant
        # (Q,R)=Q_R_maker(fitnesses_data)
        (Q,R)=Q_R_maker_self_loop(fitnesses_data)
        (Q,R)=Q_R_normalizer(Q,R)
        (Phi, Tau, Cond_Tau)=prob_and_time_calc(Q,R)
        
        
        
        for node_c in range(N_nodes):
            fix_prob[config_c, node_c] = Phi[node_c, 1]
            abs_time[config_c, node_c] = Tau[node_c]
            fix_time[config_c, node_c] = Cond_Tau[node_c, 1]
            
            
        fix_prob_avg[config_c,0] = np.mean(fix_prob[config_c, :])
        abs_time_avg[config_c,0] = np.mean(abs_time[config_c, :])
        fix_time_avg[config_c,0] = np.sum(fix_time[config_c, :] * fix_prob[config_c, :])/np.sum(fix_prob[config_c, :])
        
        np.savetxt('fix_prob.csv', fix_prob, delimiter=',')
        np.savetxt('fix_time.csv', fix_time, delimiter=',')
        np.savetxt('abs_time.csv', abs_time, delimiter=',')
        
        np.savetxt('fix_prob_avg.csv', fix_prob_avg, delimiter=',')
        np.savetxt('fix_time_avg.csv', fix_time_avg, delimiter=',')
        np.savetxt('abs_time_avg.csv', abs_time_avg, delimiter=',')
        
    return   

def scenario_0_func():
    
    FIX_PROB_A=np.zeros((len(sA_vec),len(sB_vec)))
    FIX_PROB_B=np.zeros((len(sA_vec),len(sB_vec)))
    ABS_TIME=  np.zeros((len(sA_vec),len(sB_vec)))
    FIX_TIME_A=np.zeros((len(sA_vec),len(sB_vec)))
    FIX_TIME_B=np.zeros((len(sA_vec),len(sB_vec)))
    EXT_TIME_A=np.zeros((len(sA_vec),len(sB_vec)))
    EXT_TIME_B=np.zeros((len(sA_vec),len(sB_vec)))
    
    period=periods[0]
    
    PHI     =np.zeros(( N_nodes*(N_nodes-1) , 2*len(sA_vec)*len(sB_vec) ))
    np.savetxt('Phi_p='+str(period)+'.txt', PHI, delimiter=',')
    TAU     =np.zeros(( N_nodes*(N_nodes-1) , len(sA_vec)*len(sB_vec) ))
    np.savetxt('Tau_p='+str(period)+'.txt', TAU, delimiter=',')
    COND_TAU=np.zeros(( N_nodes*(N_nodes-1) , 2*len(sA_vec)*len(sB_vec) ))
    np.savetxt('Cond_Tau_p='+str(period)+'.txt', COND_TAU, delimiter=',')
    
    
    for sA_c in range(len(sA_vec)):
        sA=sA_vec[sA_c]
        for sB_c in range(len(sB_vec)):
            sB=sB_vec[sB_c]
            
            fitness_A=np.zeros(N_nodes)
            fitness_B=np.zeros(N_nodes)
            ind=0
            for node_c in range(N_nodes):
                
                fitness_A[node_c]=rA+sA*(-1)**np.floor(ind/period)
                fitness_B[node_c]=rB+sB*(-1)**np.floor(ind/period)
                
                ind+=1
            
            
            fitness_B_temp=fitness_B.copy()
            for n_c in range(N_nodes):
                fitness_B[n_c]=fitness_B_temp[(n_c+phase_shift)%N_nodes]
            
            
            fitness_A = np.loadtxt
            fitnesses_data=np.vstack([fitness_B,fitness_A]) #row[0]: resident, row[1]: mutant
            (Q,R)=Q_R_maker(fitnesses_data)
            (Q,R)=Q_R_normalizer(Q,R)
            (Phi, Tau, Cond_Tau)=prob_and_time_calc(Q,R)
            
            n=(sA_c)*len(sB_vec)+(sB_c+1)
            w=n-1
            
            PHI=np.loadtxt('Phi_p='+str(period)+'.txt', delimiter=',')
            PHI[:,2*w:2*w+2]=Phi
            np.savetxt('Phi_p='+str(period)+'.txt', PHI, delimiter=',')
            
            TAU=np.loadtxt('Tau_p='+str(period)+'.txt', delimiter=',')
            # TAU[:,s_c]=np.reshape(Tau, (  N_nodes*(N_nodes-1)  ,1) )
            TAU[:,w]=Tau
            np.savetxt('Tau_p='+str(period)+'.txt', TAU, delimiter=',')
            
            COND_TAU=np.loadtxt('Cond_Tau_p='+str(period)+'.txt', delimiter=',')
            COND_TAU[:,2*w:2*w+2]=Cond_Tau
            np.savetxt('Cond_Tau_p='+str(period)+'.txt', COND_TAU, delimiter=',')
            
            print('sA_c: '+str(sA_c)+', '+'sB_c: '+str(sB_c))
    
    return

def built_in_pp():
    
    def intersect_finder(x_plot, y_plot, value):
        
        if y_plot[0]<value+(1e-8) and y_plot[-1]<value+(1e-8):
            pass
        else:
            return
        index =0
        
        while 1:
            if index==len(y_plot)-1:
                return
            elif (y_plot[index]-value)*(y_plot[index+1]-value)<-(1e-10):
                break
            else:
                index+=1
            
        m = (x_plot[index+1]-x_plot[index])/(y_plot[index+1]-y_plot[index])
        x_star = x_plot[index]+m*(1-y_plot[index])
        
        return x_star
    
    rA_vec = np.loadtxt('rA_vec.csv', delimiter=',')
    rA_vec_tot = np.loadtxt('rA_vec_tot.csv', delimiter=',')
    mA_vec_tot = np.loadtxt('mA_vec_tot.csv', delimiter=',')
    
    fix_prob = np.loadtxt('fix_prob.csv', delimiter=',')
    fix_time = np.loadtxt('fix_time.csv', delimiter=',')
    abs_time = np.loadtxt('abs_time.csv', delimiter=',')
    
    fix_prob_avg = np.loadtxt('fix_prob_avg.csv', delimiter=',')
    fix_time_avg = np.loadtxt('fix_time_avg.csv', delimiter=',')
    abs_time_avg = np.loadtxt('abs_time_avg.csv', delimiter=',')
    
    # selected_r_set = [0.9, 1.0, 1.1] # for all
    
    # selected_inices = [0, 1, 2] #for N=3
    # selected_inices = [0, 1, 2, 3, 4] #for N=5
    # selected_inices = [0, 1, 3, 5, 6] #for N=7
    # selected_inices = [0, 1, 3, 4, 6, 7] #for N=8
    # selected_inices = [0, 1, 4, 7, 8] #for N=9
    # selected_inices = [0, 2, 8, 14, 16] #for N=17
    # selected_inices = [0, 2, 16, 17, 29, 31] #for N=32
    # selected_inices = [0, 4, 16, 28, 32] #for N=33
    # selected_inices = [0, 8, 32, 56, 64] #for N=65
    
    try:
        selected_inices = selected_inices_dict[N_nodes]
    except:
        selected_inices = [0]
        selected_inices.append(int(round(0.25*N_nodes)))
        selected_inices.append(int(round(0.50*N_nodes)))
        selected_inices.append(int(round(0.75*N_nodes)))
    selected_inices.append(N_nodes-1)
    
    fp_plot_list = []
    ft_plot_list = []
    at_plot_list = []
    
    line_c = 0
    line_c_mat = np.zeros((len(rA_vec),2), dtype=int)
    for rA_c in range(len(rA_vec)):
        rA = rA_vec[rA_c]
        line_c_mat[rA_c, 0] = line_c
        while (np.abs(rA-rA_vec_tot[line_c])<1e-6):
            line_c+=1
            if line_c==len(rA_vec_tot):
                break
        line_c_mat[rA_c, 1] = line_c
        
    np.savetxt('line_c_mat.csv', line_c_mat, fmt='%d', delimiter=',')
    
    x_star_list=[]
    rA_x_star_list=[]
    
    for rA_c in range(len(rA_vec)):
        
        rA = rA_vec[rA_c]
        fix_prob_rA = fix_prob[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1], :]
        fix_time_rA = fix_time[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1], :]
        abs_time_rA = abs_time[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1], :]
        
        if np.abs(rA-1)<1e-8:
            base_line = 1/N_nodes
        else:
            base_line = (1-1/rA)/(1-(1/rA)**N_nodes)
        
        data = fix_prob_rA/base_line
        # plotting the heatmap
        yticklabels_list = []
        for i in mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]:
            yticklabels_list.append(str(i))
        hm = sns.heatmap(data=data,
                        annot=False,  yticklabels=yticklabels_list)
        hm.invert_yaxis()
        plt.xlabel('node index (i)')
        plt.ylabel(r'$m_A$')
        plt.yticks(fontsize=7)
        # plt.ylim(4,0)
        title_text = 'normalized local fixation probability '+'\n'+ r'$(\rho_A(i, m_A, r_A)/\rho_{Moran}(r_A))$' + '\n' + 'on a line with end self loops, for ' +r'$r_A = $'+str(rA)
        plt.title(title_text)
        plt.tight_layout()
        file_name = 'FP_local_rA='+format(round(rA, 2), '.2f')+'.png'
        plt.savefig(file_name, dpi=300)
        plt.show()
        
        if rA in selected_r_set:
            plt.figure()
            for selected_ind_c in range(len(selected_inices)):
                
                selected_ind = selected_inices[selected_ind_c]
                y_plot = fix_prob_rA[:,selected_ind]/base_line
                x_plot = mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
                plt.plot(x_plot, y_plot, 'o--', label='i ='+str(int(selected_ind)), markersize=1, linewidth=0.6)
                
            plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
            plt.grid()
            title_text = 'local normalized fixation probability '+ r'$(\rho_A(m_A, r_A, i)/\rho_{Moran}(r_A))$' + '\n on a line with end self loops' +'\n '+r'$r_A = $'+str(rA)
            plt.title(title_text)
            plt.xlabel(r'$m_A$')
            plt.ylabel(r'$\rho_A(m_A, r_A, i)/\rho_{Moran}(r_A)$')
            plt.tight_layout()
            file_name = 'FP_local_rA='+format(round(rA, 2), '.2f')+'_second.png'
            plt.savefig(file_name, dpi=300)
            plt.show()
                
        
        data = np.log10(fix_time_rA)
        # plotting the heatmap
        yticklabels_list = []
        for i in mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]:
            yticklabels_list.append(str(i))
        hm = sns.heatmap(data=data,
                        annot=False,  yticklabels=yticklabels_list)
        hm.invert_yaxis()
        plt.xlabel('node index (i)')
        plt.ylabel(r'$m_A$')
        plt.yticks(fontsize=7)
        title_text = 'local fixation time '+'\n'+ r'$(log10 (\tau_A(i, m_A, r_A)) )$' + '\n' + 'on a line with end self loops, for ' +r'$r_A = $'+str(rA)
        plt.title(title_text)
        plt.tight_layout()
        file_name = 'FT_local_rA='+format(round(rA, 2), '.2f')+'.png'
        plt.savefig(file_name, dpi=300)
        plt.show()
        
        
        fix_prob_avg_rA = fix_prob_avg[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        fix_time_avg_rA = fix_time_avg[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        abs_time_avg_rA = abs_time_avg[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        
        fp_plot_list_element = np.zeros((2, len(fix_prob_avg_rA)))
        ft_plot_list_element = np.zeros((2, len(fix_time_avg_rA)))
        at_plot_list_element = np.zeros((2, len(abs_time_avg_rA)))
        
        fp_plot_list_element[0,:] = mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        ft_plot_list_element[0,:] = mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        at_plot_list_element[0,:] = mA_vec_tot[line_c_mat[rA_c, 0]:line_c_mat[rA_c, 1]]
        
        fp_plot_list_element[1,:] = fix_prob_avg_rA
        ft_plot_list_element[1,:] = fix_time_avg_rA
        at_plot_list_element[1,:] = abs_time_avg_rA
        
        fp_plot_list.append(fp_plot_list_element)
        ft_plot_list.append(ft_plot_list_element)
        at_plot_list.append(at_plot_list_element)
        
      
    plt.figure()
    for rA_c in range(len(rA_vec)):
        rA = rA_vec[rA_c]
        
        if np.abs(rA-1)<1e-6:
            base_line = 1/N_nodes
        else:
            base_line = (1-1/rA)/(1-(1/rA)**N_nodes)
            
        x_plot = fp_plot_list[rA_c][0]
        y_plot = fp_plot_list[rA_c][1]/base_line
        
        x_star = intersect_finder(x_plot, y_plot, 1)
        x_star_list.append(x_star)
        rA_x_star_list.append(rA)
        plt.plot(x_plot, y_plot, 'o--', label=r'$r_A = $'+str(rA), markersize=2, linewidth=1)
    
    plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
    plt.grid()
    title_text = 'normalized fixation probability '+ r'$(\rho_A(m_A, r_A)/\rho_{Moran}(r_A))$' + '\n on a line with end self loops'
    plt.title(title_text)
    plt.xlabel(r'$m_A$')
    plt.ylabel(r'$\rho_A(m_A, r_A)/\rho_{Moran}(r_A)$')
    plt.tight_layout()
    file_name = 'FP_vs_mA.png'
    plt.savefig(file_name, dpi=300)
    plt.show()
    
    plt.figure()
    x_plot = rA_x_star_list
    y_plot = x_star_list
    plt.scatter(x_plot, y_plot)
    # plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
    plt.grid()
    title_text = r'$m (\rho^*=1)$'
    plt.title(title_text)
    plt.xlabel(r'$r_A$')
    plt.ylabel(r'$m (\rho^*=1)$')
    plt.tight_layout()
    file_name = 'm_rho_1.png'
    plt.savefig(file_name, dpi=300)
    plt.show()
    
    plt.figure()
    for rA_c in range(len(rA_vec)):
        rA = rA_vec[rA_c]
            
        x_plot = ft_plot_list[rA_c][0]
        y_plot = np.log10(ft_plot_list[rA_c][1])
        plt.plot(x_plot, y_plot, 'o--', label=r'$r_A = $'+str(rA), markersize=2, linewidth=1)
    
    plt.legend(loc='center right', bbox_to_anchor=(1.3, 0.5))
    plt.grid()
    title_text = 'fixation time '+ r'$(log10 (\tau_A(m_A, r_A)) )$' + '\n on a line with end self loops'
    plt.title(title_text)
    plt.xlabel(r'$m_A$')
    plt.ylabel(r'$log10(\tau_A(m_A, r_A))$')
    plt.tight_layout()
    file_name = 'FT_vs_mA.png'
    plt.savefig(file_name, dpi=300)
    plt.show()

    
    return

fitness_data = np.loadtxt('fitness_data.csv', delimiter=',')
N_nodes = np.shape(fitness_data)[1]
# N_nodes = 32
np.savetxt('N_nodes.txt', np.array([N_nodes]), fmt='%d')
N_configs = np.shape(fitness_data)[0]
np.savetxt('N_configs.txt', np.array([N_configs]), fmt='%d')


# rA=1
# rB=1
# np.savetxt('rArB.txt', np.array([rA, rB]))
# r_min=min(rA,rB)

scenario_no=1
# np.savetxt('scenario_no.txt', np.array([scenario_no]), fmt='%d')

# phase_shift=0
# np.savetxt('phase_shift.txt', np.array([phase_shift]), delimiter=',', fmt='%d')

if scenario_no==0:
    delta_s=0.05
    margin=0.1
    sA_vec=np.arange(-(r_min-margin),r_min-margin+delta_s, delta_s)
    # sA_vec=np.array([0.05,0.06,0.07])
    sB_vec=sA_vec.copy()
    # sB_vec=np.array([0, 0.1])
    
    np.savetxt('sA_vec.txt', sA_vec, delimiter=',')
    np.savetxt('sB_vec.txt', sB_vec, delimiter=',')
    kxjvnkxvjnkxndkjbdkjbjdbjdbjdbjdfbvjdfbjhbhj
    periods=[4]
    np.savetxt('periods.txt', periods, delimiter=',', fmt='%d')
    
    scenario_0_func()
    
else:
    
    # s_val=np.array([ 0, 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95])
    # s_val=np.arange(0.8,0.99,0.01)
    # s_val=np.arange(0,min(rA,rB)-0.1, 0.025)
    # s_val=np.arange(0,0.8, 0.025)
    # delta_s=0.05
    # margin=0.1
    # s_val=np.arange(0,0.8, 0.025)
    # s_val=np.arange(0,np.float32(r_min-margin+delta_s), delta_s)
    # s_val=np.arange(0,np.float32(r_min-margin+delta_s)-(1e-6), delta_s)
    # s_val=np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    # s_val=np.array([0.2, 0.4, 0.6, 0.8])
    # s_val=np.arange(0,0.4, delta_s)
    # s_val=[0, 0.19, 0.38, 0.57, 0.76, 0.95]
    # np.savetxt('s_val.txt', s_val, delimiter=',')
    # rA_val=np.array([0.8, 0.9, 1.0, 1.1])
    # np.savetxt('rA_val.txt', rA_val, delimiter=',')
    
    # periods=[1,3]
    # periods=[1,2,4]
    # periods=[1,2,4,8]
    # periods=[1,2,4,8,16]
    # periods=[1,2,4,8,16,32]
    # periods=[1,2,4,5,8,10,20,40]
    # periods=[1,2,5,10,25,50]
    # periods=[16,32]
    # periods=[1,2,4,8] #this is actually half perids (T/2)
    # np.savetxt('periods.txt', periods, delimiter=',', fmt='%d')
    
    
    configs=list(range(0,N_configs))
    np.savetxt('configs.txt', configs, delimiter=',', fmt='%d')
    
    selected_r_set = [0.9, 1.0, 1.1] # for all    
    
    selected_inices_dict=dict()
    selected_inices_dict[3] = [0, 1, 2] #for N=3
    selected_inices_dict[5] = [0, 1, 2, 3, 4] #for N=5
    selected_inices_dict[7] = [0, 1, 3, 5, 6] #for N=7
    selected_inices_dict[8] = [0, 1, 3, 4, 6, 7] #for N=8
    selected_inices_dict[9] = [0, 1, 4, 7, 8] #for N=9
    selected_inices_dict[17] = [0, 2, 8, 14, 16] #for N=17
    selected_inices_dict[32] = [0, 2, 16, 17, 29, 31] #for N=32
    selected_inices_dict[33] = [0, 4, 16, 28, 32] #for N=33
    selected_inices_dict[65] = [0, 8, 32, 56, 64] #for N=65
    
    
    
    scenario_non_0_func()
    built_in_pp()

