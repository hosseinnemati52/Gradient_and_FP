This is solving the problem by transition matrix method.
1d linear graph, with self-loops on the ends. Number of nodes determined in the file N_nodes.txt
Fitness of resident is 1.0 everywhere. 
The file rA_vec.csv includes the values of the inherent fitness of mutant (rA) for which you want to solve the problem.
First, run fitness_gen.py to produce fitness files and all configurations.
Then, run self_loop_heter_custom_fitness_TM.py to solve the problem using transition matrix.
The file fitness_data.csv includes the fitness of the mutant in each configuration and on each node. Each row is a configuration, and each column shows a node of the 1-d graph.
The file fix_prob.csv includes the fixation probability of the mutant, starting from each node. Each row belongs to a configuration, and each column shows a node of the 1-d graph.
The file fix_prob_avg.csv shows the average fixation probability of the mutant in each of the fitness configurations.
The quantity m in the plots, shows the gradient of fitness for the mutant. m=0 means no fitness gradient.