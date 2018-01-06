%% Graph script to plot results of all files

%% Parameters
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=320;
k2NL=2e4;

w1=sqrt((k1)/m1)/(2*pi);
w2=sqrt((k2)/m2)/(2*pi);
%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\DuffingValDec17\DuffOneUnitTrans135.csv';
M=csvread(file,1,0); %start reading from row 1, column 1

ansys_time = M((1:length(M)),1); % time
t=ansys_time;
ansys_amp_1 = M((1:length(M)),2);
% ansys_amp_2 = M((1:length(M)),3);
% velo1=M((1:length(M)),4);
% velo2=M((1:length(M)),5);

bandpassFile='C:\ANSYS\Temp\Validation\DuffingValDec17\HigherAmp.csv';
bandpass1=csvread(bandpassFile);
bandpass=bandpass1(:,2);
% [q,t5]=max(bandpass);
% figure
% plot(t,bandpass)
dt=abs(mean(diff(t))); 
%% TF estimate
% takes in input and output signal
% bandpass filter input signal

% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(bandpass,ansys_amp_1,[],[],[],1/dt);

% plot
% figure
hold on
graph=plot(frequencies/w1,20*log10(abs(txy)),'k');
% axis([0 12 -Inf Inf])
grid on
title('Transfer function of metamaterial configurations','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% axis([0 10 -Inf 0])
legend({'linear AMM','noise insulation panel','nonlinear case 1','nonlinear case 2','nonlinear case 3','nonlinear case 4'},'FontSize',14)
set(graph,'LineWidth',1.5);
% %%
% figure
% % hold on
% graph1=semilogx(frequencies/w1,abs(txy),'b');
% % axis([0 3 0 Inf])
% grid on
% % title('TF Estimate 5unitAMM NLH chain','FontSize',14)
% xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
% ylabel('Magnitude, dB','FontSize',14)
% set(gca,'fontsize',20)
% legend({'linear AMM numerical result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
% set(graph1,'LineWidth',2);
% %%
% figure
% % hold on
% [pxx,f] = periodogram(ansys_amp_1,[],[],1/dt);
% plot(f/w1,10*log10(pxx),'k')
% grid on
% title('Periodogram PSD of 10 unit cell linear system','FontSize',20)
% xlabel('Frequency, \omega','FontSize',20)
% ylabel('Magnitude, dB','FontSize',20)