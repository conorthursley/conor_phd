%% Dispersion Iteration
%
%------------------------------------------------------------------------------
% Purpose: Iteration of the spring and mass ratio to generate dispersion curves
%------------------------------------------------------------------------------
% Based on this dispersion relation from Huang et. al. 2009 "On the negative
% effective mass density in acoustic metamaterials" eq (7)
% -----------------------------
% m*n*w^4-((m+n)*l+2*n*k(1-cosQ))*w^2+2*k*l*(1-cosQ)=0
% --------------------------------------
% where m=m1; n=m2; k=k1; l=k2; Q=qL (wavenumber); w = frequency
tic

k1=2e6;
m1=1;
qL=linspace(0,pi,10); %evenly spaced vector from 0 to pi
A=qL;
kr = linspace(0,2,10);
hold on

parfor ii=1:length(kr)
    S_r=kr(ii); % spring ratio
%     w=linspace(0,25,1000);
    k2=k1*0.1;
    m2=m1*S_r;
    %w is angular freqeuncy and changes
    % w0 is our resonance frequency and is w0=sqrt(k^2/m^2)
    w0=sqrt(k2/m2); % rad/s
    % qL is the dimensionless wave number. in this case, the values will take
    % creat w12 which means roots 1 and 2 of angular frequency, w.
    % this value takes on a + and - value, hence 1 and 2
    % this value was created using Wolfram Alpha computational solver for
    % quartics with the equation used as f in the above section.
    w12=sqrt(-sqrt((-2*m1*k1*(1-cos(A))-k2*m1-k2*m2).^2-4*m1*m2*(2*k1*k2-2*k1*k2*cos(A)))/(m1*m2)+(2*k1*(1-cos(A)))/m1+k2/m1+k2/m2)/sqrt(2);
    % need to make the frequency ratio for the dispersion curve
    Omega1=abs(w12)/w0;

    % creat w34 which means roots 3 and 4 of angular frequency, w.
    % this value takes on a + and - value, hence 3 and 4
    w34=sqrt(sqrt((-2*m1*k1*(1-cos(A))-k2*m1-k2*m2).^2-4*m1*m2*(2*k1*k2-2*k1*k2*cos(A)))/(m1*m2)+(2*k1*(1-cos(A)))/m1+k2/m1+k2/m2)/sqrt(2);
    % need to make the frequency ratio for the dispersion curve
    Omega2=abs(w34)/w0;
    
    cell1(:,1,ii)=Omega1';
    cell2(:,1,ii)=Omega2';
    
end
toc

%% graph the dispersion curves
figure
hold on
for ii=1:length(kr)
    temp1 = cell1(:,1,ii);
    temp2 = cell2(:,1,ii);
    plot(A,temp1,'r',A,temp2,'b')
    
end
    
xlabel('Wave Number')
ylabel('Frequency ratio, w_o/w')
hold off
grid
toc


% 
% 
% figure
% plot(A,+Omega1,A,+Omega2)
% grid on
% xlabel('wavenumber, qL') % x-axis label
% ylabel('frequency ratio, w/w_0') % y-axis label
% title('Analytical determination of Dispersion Curve')
