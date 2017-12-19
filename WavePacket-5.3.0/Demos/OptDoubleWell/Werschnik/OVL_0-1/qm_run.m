%% Solve TISE
qm_setup;  
qm_init; 
qm_bound;            
qm_cleanup;

%% Matrix representations
qm_setup; 
qm_init; 
qm_matrix (pwd,'bound'); 
qm_cleanup;

%% Switch to ABNCD control language
qm_setup;  
qm_init; qm_abncd ('tdse');      
qm_cleanup;

%% Optimal control
qm_setup;  
qm_init; 
qm_optimal ('tdse'); 
qm_cleanup;