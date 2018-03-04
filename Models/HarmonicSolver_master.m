clear all
close all


n=1;
%---------------------------------------------------
m1=0.1;
zeta=0.00002;


k1=1000;

c1=zeta*2*sqrt(k1*m1);  
w1=sqrt(k1/m1);
%---------------------------------------------------
K=[k1];

C=[c1];
%---------------------------------------------------
% Mass matrix
M=[m1];

%---------------------------------------------------
% Eigenvalues
[V,D] = eig(K,M);
%---------------------------------------------------
% Forcing Frequency 
B=[1];
% S=signal(10000);

%---------------------------------------------------
% Frequency 
w=logspace(-1,2.5,10000);
%---------------------------------------------------
% Dynamic Stiffness matrix
% K_d=K-(w.^2)M;
% for loop to iterate around omega, w

w_norm=w;
disp=zeros(n);

for i=1:length(w_norm)
    K_d=K-(w(i)^2)*M-w(i);
    B=[1];
    %     displacements
    disp(:,i)=inv(K_d)*B';
      
end


%---------------------------------------------------
% %% Displacement Time responses
% figure
% plot(w_norm,(disp(1,:)),'b-')
% % axis([0 5 -20 20])
% grid on
% str = sprintf('Natural scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d', n,m1,k1);
% title(str,'FontSize',14)
% xlabel('Excitation frequency, \Omega','FontSize',14)
% ylabel('Magnitude','FontSize',14)
% legend({'mass_1'},'FontSize',14)
% 
%---------------------------------------------------
% %% Loglog plot of frequency response
% figure
% loglog(w,abs(disp(1,:)),'b-')
% y1=get(gca,'ylim');
% hold on
% % for ii=1:length(D)
% %   SP=D(ii,ii); %your point goes here
% %   line([SP SP],[y1],'Color','r', 'Linestyle','-')
% % end
% % for i=1:3
% %     SP=[0, w2*sqrt(1+mr), w2]; %your point goes here
% %     line([SP(i) SP(i)],[y1],'Color','r', 'Linestyle','-')
% % end
% hold off
% % axis([0.1 10 0.0001 10^4])
% grid on
% str = sprintf('LogLog scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d', n,m1,k1);
% title(str,'FontSize',14)
% xlabel('Frequency, rad/sec','FontSize',14)
% ylabel('magnitude','FontSize',14)
% legend({'mass_1'},'FontSize',14)
%---------------------------------------------------
%% SemiLogY Plot
figure
semilogy(w,abs(disp(1,:)),'b-')
grid on
% axis([0 3.5 10e-20 10e-5])
str = sprintf('Semilog scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d', n,m1,k1);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1'},'FontSize',14)

% function wn=signal(maxt)
% deltat=1/1000;  % [s] time increment for analysis
% % maxt=10000;     % maximum number of analysis time steps maxt*deltat
% % *DIM,excitef,ARRAY,maxt
% % *VFILL,excitef,RAND,-1.0,1.0  ! fill the array with random numbers between -1 to 1
%  
% % Define time vector
% t=deltat:deltat:deltat*maxt;
% t=t.';
%  
% % Generate white noise between -1 to +1
% whitenoise = -1 + (2).*rand(length(t),1);
%  
% % plot(t,whitenoise)
%  
% % used sptool to filter the data
% load filtered_white_data.txt
%  
% fs_white=1000;
% FFTsize=1024;
% [P_LP50Hz_force,Freq_LP50Hz_force]=pwelch(filtered_white_data,hanning(FFTsize),[],FFTsize,fs_white);
%  
% % figure
% % plot(Freq_LP50Hz_force,10*log10(abs(P_LP50Hz_force)));
% % xlabel('Frequency (Hz)');
% % ylabel('Amplitude (dB)');
% wn=Freq_LP50Hz_force;
% end
