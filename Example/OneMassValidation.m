clear all
close all


n=1;
%---------------------------------------------------
m1=0.1;
k1=3;
w1=sqrt(k1/m1);
%---------------------------------------------------
K=[k1];
%---------------------------------------------------
% Mass matrix
M=[m1];
%---------------------------------------------------
% Frequency 
w=logspace(-1,1,4000);
%---------------------------------------------------
% Eigenvalues
[V,D] = eig(K,M);
%---------------------------------------------------
% Forcing Frequency 
B=[1];

%---------------------------------------------------

% Dynamic Stiffness matrix
% K_d=K-(w.^2)M;
% for loop to iterate around omega, w

w_norm=w/w1;
disp=zeros(n);

for i=1:length(w_norm)
    K_d=K-(w(i)^2)*M;
    %     displacements
    disp(:,i)=inv(K_d)*B';
      
end


%---------------------------------------------------
%% Displacement Time responses
figure
plot(w_norm,abs(disp(1,:)),'b-')
% axis([0 5 -20 20])
grid on
str = sprintf('Natural scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d', n,m1,k1);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1'},'FontSize',14)

%---------------------------------------------------
%% Loglog plot of frequency response
figure
loglog(w,abs(disp(1,:)),'b-')
y1=get(gca,'ylim');
hold on
% for ii=1:length(D)
%   SP=D(ii,ii); %your point goes here
%   line([SP SP],[y1],'Color','r', 'Linestyle','-')
% end
% for i=1:3
%     SP=[0, w2*sqrt(1+mr), w2]; %your point goes here
%     line([SP(i) SP(i)],[y1],'Color','r', 'Linestyle','-')
% end
hold off
% axis([0.1 10 0.0001 10^4])
grid on
str = sprintf('LogLog scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d', n,m1,k1);
title(str,'FontSize',14)
xlabel('Frequency, rad/sec','FontSize',14)
ylabel('magnitude','FontSize',14)
legend({'mass_1'},'FontSize',14)%---------------------------------------------------
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