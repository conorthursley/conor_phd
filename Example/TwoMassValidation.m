clear all
close all


n=1;
%---------------------------------------------------
m1=1;
m2=0.3*m1;%kg 
k1=4.8*10^9;
k2=0.1*k1;  % N/m
mr=m2/m1;
kr=k2/k1;
w2=sqrt(k1/m2);
w1=sqrt(k2/m1);
%---------------------------------------------------
% first row, uses the first cell's secondary mass's displacement, y(2)
%calculate the stiffness using the piecewise function
% % populate
% K=zeros(n*2);
% K(1,1)=k2+(k1);
% K(1,2)=-k2;
% %K(1,3)=-k1;
% %---------------------------------------------------
% % last 2x2 diagonal, uses the last cell's secondary mass, y(end)
% % populate
% q=2*n;
% K(q,q)=k2;
% K(q-1,q-1)=k2+k1;
% K(q-1,q)=-k2;
% K(q,q-1)=-k2;
% K(q-1,q-2)=-k1;
% K(q-2,q-1)=-k1;
%---------------------------------------------------
% Iteration of remainder of cell, ignoring the first row and last diagonal
% for i=2:((2*n)-2)
%     if mod(i,2)==1  % odd number and therefore a primary mass
%         K(i,i)=k2+(k1);
%         K(i,i-1)=-k1;
%         K(i-1,i)=-k1;
%     elseif mod(i,2)==0 % even number and therefore a secondary mass
%         K(i,i)=k2+k1;
%         K(i,i-1)=-k2;
%         K(i-1,i)=-k2
%     end
% end
K=[9 -3;-3 3];
%---------------------------------------------------
% Mass matrix
m_vector=[m1 m2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
%---------------------------------------------------
% Frequency 
w=logspace(-1,1,4000);
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

w_norm=w;
disp=zeros(2*n);

for i=1:length(w_norm)
%     B((2*n)-1)=1;
    B=[0 1];
    K_d=K-(w(i)^2)*M;
%     displacements
    disp(:,i)=inv(K_d)*B';
      
end


%---------------------------------------------------
%% Displacement Time responses
figure
plot(w_norm,(disp(1,:)),'b-',w_norm,(disp(2,:)),'r-')
axis([0 5 -20 20])
grid on
str = sprintf('Natural scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n*2,m1,k1,m2,k2);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

%---------------------------------------------------
%% Loglog plot of frequency response
figure
loglog(w,abs(disp(1,:)),'b-',w_norm,abs(disp(2,:)),'r-')
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
str = sprintf('Log plot frequency response magnitudes for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n*2,m1,k1,m2,k2);
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
semilogy(w,abs(disp(1,:)),'b-',w,abs(disp(2,:)),'r-')
grid on
% axis([0 3.5 10e-20 10e-5])
str = sprintf('Semilog plot frequency response magnitudes for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n*2,m1,k1,m2,k2);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

