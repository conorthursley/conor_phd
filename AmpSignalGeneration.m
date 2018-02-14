%% Code to generate a csv file that contains a signal created by this script for use in APDL input excitation force
clear all

% input parameters
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=1.5*k1;


Fsteps=0.25; %Hz 0.25 Hz each up step
il=30; %interval length
A=1; %amplitude 
low=10; %low interval
high=40; %high interval
%
hf=high; %highest frequency
dt=1/(20*hf); %time step size
lt=((high-low)/Fsteps)*il;
% Fsteps=(high-low)/(lt/il); %freq step size
freq=34;
Tsteps=dt;
time=dt:Tsteps:il;
excel=[];
A=0.09:0.01:0.25;
tim=[];
%% create vector, then write to csv file.
for i=1:17 %how many intervals
        data=A(i)*sin(time.*freq*2*pi);
        
        excel=[excel data];
end
lt=i*il;
timeTotal=dt:Tsteps:lt;

% excel=fliplr(excel);
var=[timeTotal' excel'];

figure
plot(timeTotal,(excel))

file='U:\_PhD\APDL\Validation\DuffingValDec17\LinearSweepTS6.csv';
csvwrite(file,var);