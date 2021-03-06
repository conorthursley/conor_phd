clear all



n=1;
%---------------------------------------------------
m1=0.1;
m2=0.5*m1;%kg 
k1=1000;
k2=1.5*k1;  % N/m
mr=m2/m1;
kr=k2/k1;
w2=sqrt(k2/m2)/(2*pi);
w1=sqrt(k1/m1)/(2*pi);
zeta1=0.0000002;
zeta2=0.0000002;
c1=zeta1*2*sqrt(k1/m1);
c2=zeta2*2*sqrt(k2/m2);

%---------------------------------------------------
% damping matrix
C=zeros(n*2);
C(1,1)=c2+2*c1;
C(1,2)=-c2;
C(2,1)=-c2;
C(2,2)=c2;

%---------------------------------------------------
% first row, uses the first cell's secondary mass's displacement, y(2)
%calculate the stiffness using the piecewise function
% % populate
K=zeros(n*2);
K(1,1)=k2+(2*k1);
K(1,2)=-k2;
% If n is 1
K(2,1)=-k2;
K(2,2)=k2;


if n>1
    K(1,3)=-k1;
    K(3,1)=-k1;
else
     
end

% %---------------------------------------------------
% last 2x2 diagonal, uses the last cell's secondary mass, y(end)
% populate
if n>1
    q=2*n;
    K(q,q)=k2;
    K(q-1,q-1)=2*k1+k2;
    K(q-1,q)=-k2;
    K(q,q-1)=-k2;

end

%---------------------------------------------------
% Iteration of remainder of cell, ignoring the first row and last diagonal
for i=2:((2*n)-2)
    if mod(i,2)==1  % odd number and therefore a primary mass
        K(i,i)=k2+(2*k1);
        K(i,i+2)=-k1;
        K(i+2,i)=-k1;
    elseif mod(i,2)==0 % even number and therefore a secondary mass
        K(i,i)=k2;
        K(i,i-1)=-k2;
        K(i-1,i)=-k2;
    end
end
% K=[9 -3;-3 3];
%---------------------------------------------------
% Mass matrix
m_vector=[m1 m2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
%---------------------------------------------------
% Frequency 
w=logspace(-1,3,10000);
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
    B((2*n)-1)=0.1*k1;
%     B=[0 0.1*k1*sin(1.8)];
    K_d=K-(w(i)^2)*M-(w(i))*C;
%     displacements
    disp(:,i)=inv(K_d)*B';
      
end

x1=disp(1,:);
x2=disp(2,:);
%---------------------------------------------------
%% Displacement Time responses
figure
plot(w_norm/(w2*2*pi),x1,'b-',w_norm/(w2*2*pi),x2,'r-')
% axis([0 5 -20 20])
grid on
str = sprintf('Natural scale frequency response function for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n,m1,k1,m2,k2);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

%---------------------------------------------------
%% Loglog plot of frequency response
figure
loglog(w/(w2*2*pi),abs(x1),'b-',w_norm/(w2*2*pi),abs(x2),'r-')
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
str = sprintf('Log plot frequency response magnitudes for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n,m1,k1,m2,k2);
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
semilogy(w/(w2*2*pi),abs(disp(1,:)),'b-',w/(w2*2*pi),abs(disp(2,:)),'r-')
grid on
% axis([0 3.5 10e-20 10e-5])
str = sprintf('Semilog plot frequency response magnitudes for a %d mass-spring undamped system with m_1=%d, k_1=%d, m_2=%d, and k_2=%d', n,m1,k1,m2,k2);
title(str,'FontSize',14)
xlabel('Excitation frequency, \Omega','FontSize',14)
ylabel('Magnitude','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)

% tfestimate(x1,x2)