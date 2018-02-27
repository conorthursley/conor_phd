%% Short-time simulation
cd 1;   
qm_setup(); 
qm_init(); 
qm_propa(); 
qm_cleanup();   
cd ..;

%% Long-time simulation
cd 2;   
qm_setup(); 
qm_init(); 
qm_propa(); 
qm_cleanup();   
cd ..;
