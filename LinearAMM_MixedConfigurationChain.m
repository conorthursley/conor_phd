%% Text file writer
% creates a text file/inp file from MATLAB which is then run in ANSYS APDL
strTitle='LinearAMM_MixedConfigurationChain.txt';
fileID = fopen(strTitle,'w');
% dont forget to change the last line when you change this title
%-----------------------------------------------
%% Future Improvements and plans
% - assign variables for elements so changing between types and real
% constants is easy af.

%-----------------------------------------------
%% Intro comments and time stamp
dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
strIntro='! Script to examine the response for a linear mass-spring system\n! Conor MacDonald ';
fprintf(fileID,strcat(strIntro,dt));
%-----------------------------------------------
%% Clear up and start Proprocessor
strFIN='\n\nFINI\n/clear,all\n! Start the preprocessor\n/prep7\n';
fprintf(fileID,strFIN);
%-----------------------------------------------
%% Parameters and Modelling
% Parameters and loop to create the model
% stiffness, mass, number of cells, length, etc
%----------------10 HZ---------------------
m1=0.1429; T1=1; %element and constant type of 1
m2=0.3; T2=2; %element and constant type of 2
k1=238.69; T3=3; %element and constant type of 3
k2=118.44; T4=4; %element and constant type of 4
%----------------20 HZ---------------------
k5=961.71; T5=5; %element and constant type of 5
k6=473.74; T6=6; %element and constant type of 6
%----------------40 HZ---------------------
m7=0.2927; T7=7; %element and constant type of 7
k8=8557; T8=8; %element and constant type of 8
k9=1895; T9=9; %element and constant type of 9
%----------------70 HZ---------------------
m10=0.5176; T10=10; %element and constant type of 10
k11=47937; T11=11; %element and constant type of 11
k12=5803; T12=12; %element and constant type of 12
%-------------------------------------
L=2; %length between unit cells (cell is two masses)
l=L/2; %length within each cell
%-------------------------------------
n=100; %number of cells, so we need 2xn number of nodes
%-------------------------------------
ival=0; %initial value for node generation
fval=2*n; %final value for end of node chain
%-------------------------------------

%-------------------------------------
strPAR=('\n! Define parameters\nm1=%f\nm2=%f\nk1=%f\nk2=%f\nL=%f\nl=%f\nn=%f !number of cells\nival=%f\nfval=%f\n');
fprintf(fileID,strPAR,m1,m2,k1,k2,L,l,n,ival,fval);
%-----------------------------------------------
%% Defining elements
fprintf(fileID,'\n! --------------Elements------------------');
%----------------10 HZ---------------------
fprintf(fileID,'\n! --------------10 HZ------------------');
% Define the mass elements
fprintf(fileID,'\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,%d\n',T1,T1,m1);
fprintf(fileID,'\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,%d\n',T2,T2,m2);
% Define the linear spring elements
fprintf(fileID,'\n! Define the spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T3,T3,k1);
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T4,T4,k2);
%----------------20 HZ---------------------
fprintf(fileID,'\n! --------------20 HZ------------------');
% Define the linear spring elements
fprintf(fileID,'\n! Define the spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T5,T5,k5);
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T6,T6,k6);
%----------------40 HZ---------------------
fprintf(fileID,'\n! --------------40 HZ------------------');
% Define the mass elements
fprintf(fileID,'\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,%d\n',T7,T7,m7);
% Define the linear spring elements
fprintf(fileID,'\n! Define the spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T8,T8,k8);
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T9,T9,k9);
%----------------70 HZ---------------------
fprintf(fileID,'\n! --------------70 HZ------------------');
% Define the mass elements
fprintf(fileID,'\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,%d\n',T10,T10,m10);
% Define the linear spring elements
fprintf(fileID,'\n! Define the spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T11,T11,k11);
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,%d\n',T12,T12,k12);
%-----------------------------------------------
%-----------------------------------------------
%% Node Generation
fprintf(fileID,'\n! --------------Node Generation------------------');
% Outline the DO loop to create 'n' amount of nodes
fprintf(fileID,'\n! Define a DO loop using the command\n!*DO, counter, InitialVAL, FinalVAL, INCrement\n');
fprintf(fileID,'\n*DO,II,%d,fval,1 ! For I = %d to %d:\n',ival,ival,fval);
% Define the node command line
% we aren't concerend with rotation so we just need the longitudinal value
fprintf(fileID,'\n! Define the nodes using the command\n! N, NODE, X, Y, Z, THXY, THYZ, THZX\n');
fprintf(fileID,'\nposx=II*l		! calculate nodal position with spacing, =%d\n',l);
% define the node(s) and end the do loop
fprintf(fileID,'\nN,II+1,posx,0,0			! define the node\n*ENDDO\n');
% fprintf(fileID,'\nN,II+2,posx+l,0,0\n*ENDDO\n');
% define first node
% fprintf(fileID,'\n\nN,1,0,0,0     !Define first node at 0,0');
fprintf(fileID,'\n! Now that all the nodes are defined\n! One can define the elements that link them together\n');
%-----------------------------------------------
fprintf(fileID,'\n! --------------Element Modelling------------------');
% Masses to secondary mass
fprintf(fileID,'\n! Secondary mass is constant throughout. Do first');
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T2,T2,T2,T2);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',ival+3,fval+1);
fprintf(fileID,'\nE,II\n*ENDDO\n');
%-----------------------------------------------
fprintf(fileID,'\n! --------------10 HZ Unit Cell------------------');
% place masses on each node in the section + the next section as they are
% both the same mass (10 and 20Hz)
fprintf(fileID,'\n! Mass 1 for the 10Hz and 20Hz section are the same value, so Do loop for both');
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T1,T1,T1,T1);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',ival+2,(fval/2));
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II\n*ENDDO\n');
% link springs together (linear springs - primary mass to primary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T3,T3,T3,T3);
fprintf(fileID,'\nE,1,2\n !Define the first node to just be a mass (that will be grounded) and a primary spring to the chain');
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC ! Setout the DO loop for the first %d cells\n*DO, II,2,%d, 2\n',fval/8,(fval/4)-2);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+2\n*ENDDO\n');
%link secondary springs together (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T4,T4,T4,T4);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,2,%d, 2\n',fval/4);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
fprintf(fileID,'\n! --------------20 HZ Unit Cell------------------');
% link springs together (linear springs - primary mass to primary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T5,T5,T5,T5);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC ! Setout the DO loop for the second %d cells\n*DO, II,%d,%d, 2\n',fval/8,(fval/4),(fval/2)-2);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+2\n*ENDDO\n');
% link 20 hz springs k6, T6 (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T6,T6,T6,T6);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval/4)+2,(fval/2));
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
fprintf(fileID,'\n! --------------40 HZ Unit Cell------------------');
% primary mass
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T7,T7,T7,T7);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',(fval/2)+2,(fval*0.75));
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II\n*ENDDO\n');
% link springs together (linear springs - primary mass to primary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T8,T8,T8,T8);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC ! Setout the DO loop for the third %d cells\n*DO, II,%d,%d, 2\n',fval/8,fval/2,(fval*0.75)-2);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+2\n*ENDDO\n');
%link secondary springs together (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T9,T9,T9,T9);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval/2)+2,fval*0.75);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
fprintf(fileID,'\n! --------------70 HZ Unit Cell------------------');
% primary mass
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T10,T10,T10,T10);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',fval*0.75+2,(fval));
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II\n*ENDDO\n');
%link secondary springs together (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T12,T12,T12,T12);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval*0.75)+2,fval);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
% link springs together (linear springs - primary mass to primary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T11,T11,T11,T11);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC ! Setout the DO loop for the last %d cells\n*DO, II,%d,%d, 2\n',fval/8,fval*0.75,fval);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+2\n*ENDDO\n');
%-----------------------------------------------
%-----------------------------------------------
%% Loads: Define the excitation force on node#1
fprintf(fileID,'\n! --------------Loads/Constrain DOF------------------');
% Force of 1 applied on the first node
% prepare for harmonic or transient analysis
% fprintf(fileID,'\n! define the excitation force on Node #1, in the X direction\n! F, NODE, Lab, VALUE, VALUE2, NEND, NINC');
% fprintf(fileID,'\n! Specifies force loads at nodes.');
% fprintf(fileID,'\nF,1,FX,1');
%-----------------------------------------------
% Constrain the other end of the chain
fprintf(fileID,'\n! Constrain the first node, which is numbered %d\n! D, Node, Lab, VALUE, VALUE2, NEND, NINC, Lab2, Lab3, Lab4, Lab5, Lab6',ival+1);
fprintf(fileID,'\n! Defines degree-of-freedom constraints at nodes.');
fprintf(fileID,'\nD,%d,UX,0\n',ival+1);
%-----------------------------------------------
fprintf(fileID,'\nSAVE\nFINI\n');
%-----------------------------------------------
% Model process has been completed
%% Solution of system
% setup the solution process
% fprintf(fileID,'\n/SOLU     ! Start the solution process\n');
% fprintf(fileID,'\nANTYPE, HARMIC   ! Use Harmonic Analysis\n');
% fprintf(fileID,'!Setup the solution process\n!*\n!*\nHROPT,FULL\nHROUT,ON\nLUMPM,1\n!*\nEQSLV, ,0,\nPSTRES,0\n!*\n!*\nOUTPR,ALL,ALL,\n');
% fprintf(fileID,'HARFRQ,%d,%d,',ivalHarm,fvalHarm);
% fprintf(fileID,'\nNSUBST,%d,\nKBC,0\n!*\nSAVE\n',subsN);
% %-----------------------------------------------
% %solve the system
% fprintf(fileID,'\nSOLVE\n');
%-----------------------------------------------
%% Analysis/PostProcessing
%-----------------------------------------------
% Not working great last time. Might have to do it manually from here. 
%-----------------------------------------------
% Plot the graph of the response of first and last node in the chain
% fprintf(fileID,'\n! Plot the graph of the response of first and last node in the chain');
% fprintf(fileID,'\n/POST26\nFILE,''file'',''rst'',''.''\nNUMVAR,200\n!SOLU,191,NCMIT\n');
% fprintf(fileID,'\n!STORE,MERGE\nPLCPLX,0\nPRCPLX,1\n!*\n');
% fprintf(fileID,'\n! Define results set 2, to be the UX displacement of Node #1\n! and give it a text label UX_1 to print on the graph\n');
% fprintf(fileID,'NSOL,2,1,U,X, UX_1,\nNSOL,3,chainlen+1,U,X, UX_end,\nXVAR,1	! Define the x label to be the frequency (set 1)\n');
% fprintf(fileID,'PLVAR,2,3');


%% End
% Copy text file contents to ANSYS' input file format, .inp
fclose(fileID);
fileID2 = fopen('APDL_AMM_MixedLinearChain.inp' ,'w');
copyfile LinearAMM_MixedConfigurationChain.txt APDL_AMM_MixedLinearChain.inp



