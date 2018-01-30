%% Code to generate a csv file that contains a signal created by this script for use in APDL input excitation force
clear all

% input parameters
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=1.5*k1;


Fsteps=0.25; %Hz 0.25 Hz each up step
il=10; %interval length
A=1; %amplitude 
low=10; %low interval
high=40; %high interval
%
hf=high; %highest frequency
dt=1/(20*hf); %time step size
lt=((high-low)/Fsteps)*il;
% Fsteps=(high-low)/(lt/il); %freq step size
Tsteps=dt;
time=dt:Tsteps:il;
excel=[];
%% create vector, then write to csv file.
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

file='U:\_PhD\APDL\Validation\DuffingValDec17\LinearSweepTS3.csv';
csvwrite(file,var);
%% Sweep Plot

file = 'U:\_PhD\APDL\Validation\DuffingValDec17\DuffOneUnitTrans265.csv';
M=csvread(file,1,0); %start reading from row 3, column 1

time = M((1:end),1);      % excitation frequency
velo1=M((1:end),4);
velo2=M((1:end),5);

w0=(zeros(length(time),1))';
w1=(zeros(length(time),1))';
w2=w1;

freq=low+Fsteps:Fsteps:high;
ranger=1:(length(M)/(lt/il));
step_freq=low;
for ii=1:(lt/il)-1 %frequency intervals
    w0(:,ranger)=step_freq;
    mu=mean(velo1(ranger));
    mu2=mean(velo2(ranger));
    w1(:,ranger)=mu;   % average of velo1
    w2(:,ranger)=mu2;
    step_freq=step_freq+Fsteps;
    ranger=ranger+round((length(M)/(lt/il)));
end

KE1hat=(0.5*m1*w1.^2);
KE2hat=(0.5*m2*w2.^2);

RDRhat=KE2hat./(KE1hat+KE2hat);
% RDR_avg=mean(RDRhat);
% avg=ones(length(time),1);
% avg=avg*RDR_avg;

figure
hold on
plot(w0,(RDRhat),'b-*')
grid on
xlabel ('Input Frequency, Hz','FontSize',18)
ylabel ('Energy Distribution','FontSize',18)
legend S1


