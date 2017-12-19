%% Bound states
qm_setup; 
qm_init; 
qm_bound;
qm_cleanup;

%% Matrix elements
qm_setup; 
qm_init; 
qm_matrix (pwd,'bound');
qm_cleanup;

%% TDSE: Get ABNCD matrices
qm_setup; 
qm_init; 
qm_abncd('lvne');
qm_cleanup;

%% TDSE: Optimal Control Theory
qm_setup;
qm_init;
qm_optimal('lvne');
qm_cleanup();
