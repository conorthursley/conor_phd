%% Bare state representation
cd 0-1/BareState;   
qm_setup(); 
qm_init(); 
qm_propa(); 
qm_cleanup();   
cd ../..;

%% Floquet dressed state representation
cd 0-1/Dressed02;   
qm_setup(); 
qm_init(); 
qm_propa(); 
qm_cleanup();   
cd ../..;
