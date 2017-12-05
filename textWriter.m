%% Text file writer
% creates a text file/inp file from MATLAB which is then run in ANSYS APDL
strTitle='Name_of_my_text_file.txt';
fileID = fopen(strTitle,'w');
%-----------------------------------------------
%% Intro comments and time stamp
dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
strIntro='! Script to examine the response for a mass-spring system\n! Conor MacDonald ';
fprintf(fileID,strcat(strIntro,dt));
%-----------------------------------------------
%% Clear up and start Proprocessor
strFIN='\n\nFINI\n/clear,all\n! Start the preprocessor\n/prep7\n';
fprintf(fileID,strFIN);
%-----------------------------------------------
%% Parameters and Modelling
% Parameters and loop to create the model
% stiffness, mass, number of cells, length, etc
m1=1; 
m2=0.3;
k1=1000;
k2=100;
L=1; %length between cells (cell is two masses)
l=0.5; %length within each cell
n=5; %number of cells, so we need 2xn number of nodes
%-------------------------------------
ival=0; %initial value for node generation
fval=2*n; %final value for end of node chain
%-------------------------------------
% frequency range for harmonic analysis
ivalHarm=0; %initial freq
fvalHarm=100;  %final freq
subsN=10000; %number of substeps in analysis
%-------------------------------------
strPAR=('\n! Define parameters\nm1=%f\nm2=%f\nk1=%f\nk2=%f\nL=%f\nl=%f\nn=%f !number of cells\nival=%f\nfval=%f\n');
fprintf(fileID,strPAR,m1,m2,k1,k2,L,l,n,ival,fval);
%-----------------------------------------------
% Define the mass element, element type 1
strMAS=('\n! Define the mass element\nET,1,MASS21\nKEYOPT,1,3,4\nR,1,m1\n');
fprintf(fileID,strMAS);
%-----------------------------------------------
% Define the spring element, element type 2
%keyopts, real constants, declare element type
fprintf(fileID,'\n! Define the spring element\nET,2,COMBIN14\nKEYOPT,2,3,2\nR,2,k1\n');

%-----------------------------------------------
% Outline the DO loop to create 'n' amount of nodes
fprintf(fileID,'\n! Define a DO loop using the command\n!*DO, counter, InitialVAL, FinalVAL, INCrement\n');
fprintf(fileID,'\n*DO,II,%d,fval,1 ! For I = %d to %d:\n',ival,ival+1,fval);

% Define the node command line
% we aren't concerend with rotation so we just need the longitudinal value
fprintf(fileID,'\n! Define the nodes using the command\n! N, NODE, X, Y, Z, THXY, THYZ, THZX\n');
fprintf(fileID,'\nposx=II*L		! calculate nodal position with spacing, L=%d\n',L);
% define the node(s) and end the do loop
fprintf(fileID,'\nN,II+1,posx,0,0			! define the node\n*ENDDO\n');
fprintf(fileID,'\n! Now that all the nodes are defined\n! One can define the elements that link them together\n');
%-----------------------------------------------
% link springs together
fprintf(fileID,'\nTYPE,2! Change the element type to 2 (spring element)\nREAL,2! Change to real set 2 for the spring\n');
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,1 ,%f, 1\n',fval);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II,II+1\n*ENDDO\n');
%-----------------------------------------------
% place masses on each node=
fprintf(fileID,'\nTYPE,1! Change the element type to 1 (mass element)\nREAL,1! Change to real set 1 for the mass\n');
% setout another DO loop
fprintf(fileID,'\n!*DO, Par, IVAL, FVAL, INC\n*DO, II,%d ,%d, 1\n',ival+1,fval+1);
% element is defined by connectivity to two nodes, I and J
fprintf(fileID,'\nE,II\n*ENDDO\n');
%-----------------------------------------------
% Define the excitation force on node#1
% Force of 1 applied on the first node
% prepare for harmonic or transient analysis
fprintf(fileID,'\n! define the excitation force on Node #1, in the X direction\n! F, NODE, Lab, VALUE, VALUE2, NEND, NINC');
fprintf(fileID,'\n! Specifies force loads at nodes.');
fprintf(fileID,'\nF,1,FX,1');
%-----------------------------------------------
% Constrain the other end of the chain
fprintf(fileID,'\n! Constrain the end node, which is numbered %d\n! D, Node, Lab, VALUE, VALUE2, NEND, NINC, Lab2, Lab3, Lab4, Lab5, Lab6',fval);
fprintf(fileID,'\n! Defines degree-of-freedom constraints at nodes.');
fprintf(fileID,'\nD,%d,UX,0\n',fval+1);
%-----------------------------------------------
fprintf(fileID,'\nSAVE\nFINI\n');
%-----------------------------------------------
% Model process has been completed
%% Solution of system
% setup the solution process
fprintf(fileID,'\n/SOLU     ! Start the solution process\n');
fprintf(fileID,'\nANTYPE, HARMIC   ! Use Harmonic Analysis\n');
fprintf(fileID,'!Setup the solution process\n!*\n!*\nHROPT,FULL\nHROUT,ON\nLUMPM,1\n!*\nEQSLV, ,0,\nPSTRES,0\n!*\n!*\nOUTPR,ALL,ALL,\n');
fprintf(fileID,'HARFRQ,%d,%d,',ivalHarm,fvalHarm);
fprintf(fileID,'\nNSUBST,%d,\nKBC,0\n!*\nSAVE\n',subsN);
%-----------------------------------------------
%solve the system
fprintf(fileID,'\nSOLVE\n');
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
fileID2 = fopen('APDL_conor.inp' ,'w');
copyfile Name_of_my_text_file.txt APDL_conor.inp



