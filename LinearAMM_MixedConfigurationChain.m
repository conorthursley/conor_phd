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
m1=0.1; 
T1=1; %element and constant type of 1
m2=0.3;
T2=2; %element and constant type of 2
k1=1000;
T3=3; %element and constant type of 3
%-------------------------------------
% 4 types of springs. mass remains the same throughout the chain
k2=1.184e3; %translates to 10Hz - k2=0.3*(10*2*pi)^2, 0.3=m2
T4=7; %element and constant type of 4
k3=4.7374e3; %20Hz
T5=6;%element and constant type of 5
k4=1.0659e4; %30Hz
T6=5; %element and constant type of 6
k5=1.8950e4; %40Hz
T7=4; %element and constant type of 7
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
strPAR1=('\nk3=%f\nk4=%f\nk5=%f\n');
fprintf(fileID,strPAR1,k3,k4,k5);
%-----------------------------------------------
% Define the mass element, element type 1
strMAS=('\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,m1\n');
fprintf(fileID,strMAS,T1,T1);
% Define mass 2, element 2
fprintf(fileID,'\n! Define the mass element\nET,%d,MASS21\nKEYOPT,1,3,4\nR,%d,m2\n',T2,T2);
%-----------------------------------------------
% Define the linear spring element, element type 2
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the primary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,k1\n',T3,T3);
%-----------------------------------------------
% Define the secondary linear spring element, spring value that corresponds to 10Hz, element type 3
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,k2\n',T4,T4);
%-----------------------------------------------
% spring value that corresponds to 20Hz, element type 
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,k2\n',T5,T5);
%-----------------------------------------------
% spring value that corresponds to 30Hz, element type 
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,k2\n',T6,T6);
%-----------------------------------------------
% spring value that corresponds to 40Hz, element type 
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the secondary spring element\nET,%d,COMBIN14\nKEYOPT,2,3,2\nR,%d,k2\n',T7,T7);
%-----------------------------------------------
%-----------------------------------------------
%% Node Generation

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
% link springs together (linear springs - primary mass to primary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T3,T3,T3,T3);
fprintf(fileID,'\nE,1,2\n !Define the first node to just be a mass (that will be grounded) and a primary spring to the chain');
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,2,%d, 2\n',fval-2);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+2\n*ENDDO\n');
%-----------------------------------------------
% Linking 4 different spring values for K2, so need to break up the DO
% loop, or repeat it 4 times.
%------------------------------------------------
% link 10 hz springs k2, T4 (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T4,T4,T4,T4);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,2,%d, 2\n',fval/4);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
% link 20 hz springs k2, T5 (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T5,T5,T5,T5);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval/4)+2,fval/2);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
% link 30 hz springs k2, T6 (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T6,T6,T6,T6);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval/2)+2,fval*0.75);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
% link 40 hz springs k2, T7 (linear springs - primary mass to secondary mass)
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (spring element)\nREAL,%d! Change to real set %d for the spring\n',T7,T7,T7,T7);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d,%d, 2\n',(fval*0.75)+2,fval);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
% place masses on each node=
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T1,T1,T1,T1);
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',ival+2,fval);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II\n*ENDDO\n');
% Change masses to secondary mass
fprintf(fileID,'\nTYPE,%d! Change the element type to %d (mass element)\nREAL,%d! Change to real set %d for the mass\n',T2,T2,T2,T2);
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 2\n',ival+3,fval+1);
fprintf(fileID,'\nE,II\n*ENDDO\n');
%-----------------------------------------------
%% Loads: Define the excitation force on node#1

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



