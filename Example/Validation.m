    clear all

%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\MDOF_1.csv';
M=csvread(file,2,1); %start reading from row 2, column 1

ansys_freq = M((1:40000),2); 
ansys_freq_hz = M((1:40000),5); % excitation frequency
ansys_amp = M((1:40000),1); 


%% Code to setup analytical method
n=2;
%---------------------------------------------------
m1=1;
m2=0.3*m1;%kg 
k1=4.8*10^9;
k2=0.1*k1;  % N/m
mr=m2/m1;
kr=k2/k1;
w2=sqrt(k2/m2);
w1=sqrt(k1/m1);
%---------------------------------------------------
% first row, uses the first cell's secondary mass's displacement, y(2)
%calculate the stiffness using the piecewise function
% populate
K=zeros(n*2);
K(1,1)=k2+(k1);
K(1,2)=-k2;
K(1,3)=-k1;
%---------------------------------------------------
% last 2x2 diagonal, uses the last cell's secondary mass, y(end)
% populate
q=2*n;
K(q,q)=k2;
K(q-1,q-1)=k2+k1;
K(q-1,q)=-k2;
K(q,q-1)=-k2;
%---------------------------------------------------
% Iteration of remainder of cell, ignoring the first row and last diagonal
for i=2:((2*n)-2)
    if mod(i,2)==1  % odd number and therefore a primary mass
        K(i,i)=k2+(2*k1);
        K(i,i+1)=-k2;
        K(i,i+2)=-k1;
    elseif mod(i,2)==0 % even number and therefore a secondary mass
        K(i,i)=k2;
        K(i,i-1)=-k2;
        K(i+1,i-1)=-k1;
    end
end
%---------------------------------------------------
% Mass matrix
m_vector=[m1 m2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
%---------------------------------------------------
% Frequency 
w=logspace(1,8,40000); %rad/s
%---------------------------------------------------
% Eigenvalues
[V,D] = eig(K,M);
%---------------------------------------------------
% Forcing Frequency 
B=zeros(1,(n*2));

%---------------------------------------------------

% Dynamic Stiffness matrix
% K_d=K-(w.^2)M;
% for loop to iterate around omega, w

w_norm=w/w2;
disp=zeros(2*n);

for i=1:length(w_norm)
    B((2*n)-1)=1; 
    K_d=K-(w(i)^2)*M;
%     displacements
    disp(:,i)=inv(K_d)*B';
      
end


%---------------------------------------------------
%% Displacement Time responses
figure
plot(w_norm,abs(disp(1,:)),'-',w_norm,abs(disp(2,:)),'r-')
axis([0 3.5 -1e-7 1e-7])
grid on
str = sprintf('Frequency response magnitudes for a %d unit-cell AMM with m_1=%d, k_1=%d, m_2=0.3*m_1, and k_2=0.1*k_1',n, m1,k1);
title(str,'FontSize',14)
xlabel('Normalised frequency, \Omega/w_2','FontSize',14)
ylabel('Normaliesd magnitude, U_2/f','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

%---------------------------------------------------
%% Loglog plot of frequency response
figure
loglog(w_norm,abs(disp(1,:)),'-',w_norm,abs(disp(2,:)),'r-')
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
str = sprintf('Frequency response magnitudes for a %d unit-cell AMM with m_1=%d, k_1=%d, m_2=0.3*m_1, and k_2=0.1*k_1',n, m1,k1);
title(str,'FontSize',14)
xlabel('Frequency, rad/sec','FontSize',14)
ylabel('magnitude','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)
%---------------------------------------------------
% %% Normal Plot
% figure
% plot(w_norm,disp(1,:),'-')
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
% str = sprintf('Frequency response magnitudes for m_1=%d, k_1=%d, m_2=%d, and k_2=%d', m1,k1,m2,k2);
% title(str)
% xlabel('Frequency, rad/sec')
% ylabel('magnitude')
% 
% %------------------------
%% SemiLogY Plot
figure
semilogy(w_norm,abs(disp(1,:)),'-',w_norm,abs(disp(2,:)),'r-')
grid on
axis([0 3.5 10e-20 10e-5])
str = sprintf('Frequency response magnitudes for a %d unit-cell AMM with m_1=%d, k_1=%d, m_2=0.3*m_1, and k_2=0.1*k_1',n, m1,k1);
title(str,'FontSize',14)
xlabel('Normalised frequency, \Omega/w_2','FontSize',14)
ylabel('Normaliesd magnitude, U_2/f','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

%% Comparison of Numerical and Analytical
figure
loglog(w_norm,abs(disp(1,:)),'b',ansys_freq,abs(ansys_amp),'r-')
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
str = sprintf('Frequency response magnitudes for a %d unit-cell AMM with m_1=%d, k_1=%d, m_2=0.3*m_1, and k_2=0.1*k_1',n, m1,k1);
title(str,'FontSize',14)
xlabel('Frequency, rad/sec','FontSize',14)
ylabel('magnitude','FontSize',14)
legend({'matlab','APDL'},'FontSize',14)
%% TF Estimate

[txy,frequencies]=tfestimate(ansys_freq_hz,ansys_amp,[],[],[],500);

% plot
figure
plot(frequencies,20*log10(abs(txy)))
grid on
title('TF Estimate','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
