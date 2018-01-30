%% Graph script to plot results of all files

%% Parameters
m1=0.1; 
m2=0.5*m1;
k1=1000;
k2=1500;
k2NL=2e4;

w1=sqrt((k1)/m1)/(2*pi);
w2=sqrt((k2)/m2)/(2*pi);
%% File read from APDL simulation (numerical)

file = 'U:\_PhD\APDL\Validation\DuffingValDec17\DuffOneUnitTrans279.csv';
M=csvread(file,1,0); %start reading from row 1, column 1
cue=0;
ansys_time = M((1:length(M)-cue),1); % time150000:
t=ansys_time;
ansys_amp_1 = M((1:length(M)-cue),2);
ansys_amp_2 = M((1:length(M)-cue),3);
velo1=M((1:length(M)-cue),4);
velo2=M((1:length(M)-cue),5);
% if signal was generated with multiple cells and want to compare the
% difference
ff=2.393; %forcing frequency
wf1=w1*ff;
wf2=w2*ff;
% bandpass= M((1:length(M)),6);
% bandpass= 1*cos(wf1*2*pi.*t);

bandpassFile='U:\_PhD\APDL\Validation\DuffingValDec17\LinearModeTS7.csv';
bandpass1=csvread(bandpassFile);
bandpass=bandpass1(:,2);
% % bandpass=ones(1,length(M));
% bandpass=3000.*bandpass;
x1=ansys_amp_1; 
x2=ansys_amp_2;


dt=abs(mean(diff(t))); 
%% TF estimate
% takes in input and output signal
% bandpass filter input signal
wind = kaiser(length(x1),5);
% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(bandpass,x1,wind,[],[],1/dt);
figure
hold on
% subplot(2,2,4)

graph=semilogx(frequencies,20*log10(abs(txy)),'b');
% axis([0 12 -Inf Inf])
grid on
% title('Transfer function of metamaterial configurations','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% axis([0 10 -Inf 0])
% legend({'linear AMM','noise insulation panel','nonlinear case 1','nonlinear case 2','nonlinear case 3','nonlinear case 4'},'FontSize',14)
% legend({'1 unit cell','2 unit cells','5 unit cells','10 unit cells'},'FontSize',14)
set(graph,'LineWidth',1.5);
alldatacursors = findall(gcf,'type','hggroup');
set(alldatacursors,'FontSize',20)
hold off

%%
% figure
% hold on
% graph1=semilogx(frequencies/w1,abs(txy),'m');
% % axis([0 3 0 Inf])
% grid on
% % title('TF Estimate 5unitAMM NLH chain','FontSize',14)
% xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
% ylabel('Magnitude, dB','FontSize',14)
% set(gca,'fontsize',20)
% legend({'linear AMM','nonlinear case 1','nonlinear case 2','nonlinear case 3','nonlinear case 4'},'FontSize',14)
% % legend({'linear AMM numerical result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
% set(graph1,'LineWidth',2);
% alldatacursors = findall(gcf,'type','hggroup');
% set(alldatacursors,'FontSize',20)
%%
figure
% hold on
% subplot(2,2,4)

wind1 = kaiser(length(x1),15);
[pxx,f] = periodogram(x1,wind1,[],1/dt);
plot(f/w2,10*log10(pxx),'b')
grid on
% title('Periodogram PSD of 10 unit cell linear system','FontSize',20)
xlabel('Frequency, \omega/\omega_2','FontSize',20)
ylabel('Magnitude, dB','FontSize',20)
set(gca,'FontSize',24)
axis([0 2 -100 0])
