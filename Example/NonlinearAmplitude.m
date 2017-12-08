%% Validation of nonlinear mass spring system of Duffing type
% can also try van der pol oscillator later on
%% File read from APDL simulation

file = 'U:\_PhD\APDL\Duffing\FivetToSevenHz.csv';
M=csvread(file,2,19); %start reading from row 3, column 1

time = M((1:18009),1);      % excitation frequency
amp =M((1:18006),2);
w0=zeros(length(time),1);
w1=zeros(length(time),1);
%-------------------------------------------------------
% setup the freqeuncy range 5.25 to 7.25 in steps of 0.25
ranger=1:2001;
step_freq=5.25;
for ii=1:9 %8 steps from 5.25 to 7.25
    w0(ranger,1)=step_freq;
    step_freq=step_freq+0.25;
    ranger=ranger+2001;
end

    
for i=1:length(time)
    w1(i)=asin(abs(amp(i)))/(time(i))*2*pi;
    if isnan(w0(i))==1
        w1(i)=0;
    else
    end
    
end

steps=linspace(5.25,7.25,length(time));    % linearly spaced vector from f_1 to f_end
 %convert from Hz to Rad/s
%% Plots and plotting
 % Plot of U1 and U2 displacements against input frequency
figure
plot(time,amp)
grid on
xlabel 'Input Frequency, Hz'
ylabel 'Displacement, mm'
title 'Displacement of U1 and U2 as a function of w'
