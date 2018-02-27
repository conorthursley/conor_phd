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

%% LvNE: H2 optimal model reduction
dim=160;
qm_setup();
qm_init();
qm_H2model('lvne',dim); 
qm_cleanup();

%% LvNE: Optimal control theory in reduced dimensionality
qm_setup();
qm_init();
qm_optimal('lvne','h',dim);
qm_cleanup();