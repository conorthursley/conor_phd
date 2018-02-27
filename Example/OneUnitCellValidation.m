    clear all
k1=4.8e9;
m1=1;
w=sqrt(k1/m1);
%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\OneUnitCell.csv';
M=csvread(file,2,0); %start reading from row 2, column 1

ansys_freq = M((1:4000),1); % excitation frequency
ansys_amp_1 = M((1:4000),2);
ansys_amp_2 = M((1:4000),2);

%% Displacement Time responses
figure
plot(ansys_freq,ansys_amp_1)
grid on
title('Frequency response magnitudes for a one unit-cell AMM','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
% legend({'mass_1'},'FontSize',14)
%% Loglog plot of frequency response
figure
loglog(ansys_freq,(abs(ansys_amp_1)))
y1=get(gca,'ylim');
grid on
title('Frequency response magnitudes for a one unit-cell AMM','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)

