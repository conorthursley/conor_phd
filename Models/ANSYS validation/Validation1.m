%% Code to validate the MATLAB models with the ANSYS models
% 21/05/18 Conor MacDonald
tic
%% File read from APDL simulation (numerical)

file = 'U:\_PhD\ANSYS\ANSYS Validation\Completed Simulations\NLchainWithLinearParameters10seconds_TS1.csv';
M=csvread(file,1,0); %start reading from row 2, column 1
cue=0;
lenSig=1:length(M)-cue; %1:length(M)
ansys_time = M(lenSig,1); % time150000:
t=ansys_time;
ansys_amp_1 = M(lenSig-cue,2);
% ansys_amp_2 = M(lenSig-cue,3);
% ansys_amp_3 = M(lenSig,4);
% ansys_amp_4 = M(lenSig,5);

%% Displacement Time responses
figure
plot(ansys_time,(ansys_amp_1),'g.','LineWidth',0.005)
grid on
title('Time response magnitudes for a free-fixed AMM system with 2 units','FontSize',14)
xlabel('Time, s','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)
set(gca,'fontsize',14)



toc