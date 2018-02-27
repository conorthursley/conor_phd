%% Bound states
qm_setup(); 
qm_init; 
qm_bound();
qm_cleanup();

%% Matrix elements
qm_setup(); 
qm_init; 
qm_matrix ();
qm_cleanup();

%% LvNE: Get ABNCD matrices
qm_setup(); 
qm_init; 
qm_abncd('lvne');
qm_cleanup();

%% LvNE: Optimal control theory in full dimensionality
qm_setup();
qm_init();
qm_optimal('lvne');
qm_cleanup();

