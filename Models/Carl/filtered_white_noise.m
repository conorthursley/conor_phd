%% Script to generate white noise force excitation over 50 Hz

%% This is the ansys script
% ! Define an array containing random numbers
% ! to be used for the forcing function
deltat=1/1000;  % [s] time increment for analysis
maxt=10000;     % maximum number of analysis time steps maxt*deltat
% *DIM,excitef,ARRAY,maxt
% *VFILL,excitef,RAND,-1.0,1.0  ! fill the array with random numbers between -1 to 1
 
%% Define time vector
t=deltat:deltat:deltat*maxt;
t=t.';
 
%% Generate white noise between -1 to +1
whitenoise = -1 + (2).*rand(length(t),1);
 
plot(t,whitenoise)
 
% used sptool to filter the data
% load filtered_white_data.txt
filtered_white_data=whitenoise;
fs_white=1000;
FFTsize=1024;
[P_LP50Hz_force,Freq_LP50Hz_force]=pwelch(filtered_white_data,hanning(FFTsize),[],FFTsize,fs_white);
 
figure
plot(Freq_LP50Hz_force,10*log10(abs(P_LP50Hz_force)));
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
 
 
 
