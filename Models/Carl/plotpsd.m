% Plot the PSD of the results
fs = 50;    % [Hz] sampling frequency

%[Pxx,F] = pwelch(signal(:,2),hanning(length(signal)),fs);
[Pxx,F]=pwelch(signal(:,2),[],[],[],fs);
[Pxx2,F2]=pwelch(signal(:,2),hanning(length(signal)),[],[],fs);
plot(F,10*log10(Pxx),F2,10*log10(Pxx2)); xlabel('Hz'); ylabel('dB');
title('Welch''s PSD Estimate');

% Calculate the transfer function
FFTsize=1024;
force_FFT=fft(force.*hanning(length(force)),FFTsize);
disp_FFT=fft(signal(:,2).*hanning(length(force)),FFTsize);

xfer=disp_FFT./force_FFT;

df=fs/FFTsize;
freq=[0:FFTsize/2]'*df;

xfer = xfer(1:FFTsize/2+1);
xfer_mag = abs(xfer);

% Now try using the matlab 
figure

%[Txy,F3] = tfestimate(force,signal(:,2),hanning(length(signal)),[],FFTsize,fs)
[Txy,F3] = tfestimate(force,signal(:,2),1024,[],[],fs)

plot(freq,10*log10(xfer_mag),'-',F3,10*log10(abs(Txy)),'.');
xlabel('Frequency');
ylabel('FRF disp / force (dB)')
legend('CQH TF','Matlab tfestimate')





%% Matlab model of mass-spring system to compare with ansys
dt = 1/fs;
maxt = 1000;
t = dt*(0:maxt-1);

mass=1;		% [kg]
resfreq=5;	% [Hz]
stiff=mass*(2*pi*resfreq)^2;    % [N/m]
damp=0.0001;        % keep as a small number to fix solver errors

% Theoretical FRF
omega=2*pi*freq;

f_natural=sqrt(stiff/mass)/(2*pi);
ratio_freq=freq/f_natural;
zeta=damp / (2*sqrt(stiff*mass));
frf=ones(size(freq))./sqrt( (1-ratio_freq.^2).^2 + (2*zeta*ratio_freq).^2 ) / stiff;

figure
p2h=plot(freq,10*log10(abs(frf)),'r-');


Ac = [0 1; -stiff/mass -damp/mass];
Bc = [0 1/mass]';

A = expm(Ac*dt);
B = Ac\(A-eye(2))*Bc;

C = [1 0];
D = 0;


u = force;

y = 0;
x = [0;0];
for k = 1:maxt
    y(k) = C*x + D*u(k);
    x = A*x + B*u(k);
end

% figure
% plot(t,y,signal(:,1),signal(:,2),'-.');
% legend('Matlab','Ansys');
% xlabel('Time')
% ylabel('Displacement (m)')

