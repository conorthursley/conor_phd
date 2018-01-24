%% Code to generate a csv file that contains a signal created by this script for use in APDL input excitation force
clear all

% input parameters
% lt=300; %length of time
Fsteps=0.25; %Hz 0.25 Hz each up step
il=15; %interval length
A=1; %amplitude 
low=5; %low interval
high=25; %high interval
%
hf=high; %highest frequency
dt=1/(20*hf); %time step size
lt=((high-low)/Fsteps)*il;
% Fsteps=(high-low)/(lt/il); %freq step size
Tsteps=dt;
time=dt:Tsteps:il;
excel=[];
% create vector, then write to csv file.
for i=1:(lt/il) %how many intervals
        freq=low+Fsteps*(i-1);
        data=A*sin(time.*freq*2*pi);
        
        excel=[excel data];
end
timeTotal=dt:Tsteps:lt;
% excel=fliplr(excel);
var=[timeTotal' excel'];

figure
plot(timeTotal,(excel))

file='U:\_PhD\APDL\Validation\DuffingValDec17\LinearSweepTS1.csv';
csvwrite(file,var);

