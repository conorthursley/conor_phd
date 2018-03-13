%% Analytical Code for AMM model

n=7; % number of metamaterial units in the chain

%% Parameters for AMM unit cell
% parameters should be passed through ODE45 to avoid reptitive changing of
% inputs
mass1=101.1e-3;		% [kg]
mass2=46.47e-3;
stiff1=117;    % [N/m]
stiff2=74;
damp1=0.0002;     % [Ns/m] keep as a small number to fix solver errors
damp2=0.0002;
w=27*2*pi; % driving frequency

omega_axis=logspace(-1,3,4000);

%% Ratio parameters
theta=mass2/mass1;  % mass ratio 
w1=sqrt(stiff1/mass1);
w2=sqrt(stiff2/mass2);
% exctiation frequency 
w=linspace(0,200,2000);
eta_r=w/w2;       %ratio of excitation frequency to mass2 frequency 
% this value will be the independent varaible for the bandgap diagram
eta_s=w2/w1;        % structural frequency ratio

Qa=acos(1-(((1-eta_r.^2+theta)/(2*(1-eta_r.^2)))*eta_r.^2*eta_s^2));
%% Matrices of parameters
%---------------------------------------------------
% mass matrix
m_vector=[mass1 mass2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
%---------------------------------------------------
% Damping matrix
[C]=spring_matrix(n,damp1,damp2);
%---------------------------------------------------
% Stiffness matrix
[K]=spring_matrix(n,stiff1,stiff2);
%---------------------------------------------------
% Forcing Frequency 
B=[1;0]; % forcing vector of the first AMM cell

%% Matrix 
% matrix in the form seen here
% http://www.brown.edu/Departments/Engineering/Courses/En4/Notes/vibrations_mdof/vibrations_mdof.htm
N=[eye(2*n) zeros(2*n);zeros(2*n) M];
D=[zeros(2*n) -eye(2*n);K C];
f=[zeros(4*n,1)];
f(2*n+1)=1;

%% Amplitude Response loop
X=zeros(length(omega_axis),2*n)';
for jj=1:length(omega_axis)
[amp,phase] = damped_forced_vibration(D,N,f,omega_axis(jj));
X(:,jj)=amp;
end

%% Plots
% vibraiton amplitude - note: NOT TRANSMITTANCE 
figure
plot1=loglog(omega_axis,X(1,:),omega_axis,X((n-1),:));
title('Amplitude Response')
grid on
xlabel('Frequency,  (Hz)')
ylabel('Amp')
set(plot1,'LineWidth',1.5)
legend 'Outer mass' 'Inner mass'
set(gca,'fontsize',20)

% Transmission

U=X((2*n-1),:)./(X(1,:));

figure
plot2=plot(omega_axis/w2, 20*log(U));
title('Transmission')
grid on
xlabel('Frequency ratio, (Hz)')
ylabel('Transmission')

%% Dispersion Diagram
figure
plot(w/w2,real(Qa)/pi,w/w2,abs(imag(Qa)/(pi)))


%% Forced Vibraiton Amplitude function 
function [amp,phase] = damped_forced_vibration(D,M,f,omega)

% Function to calculate steady state amplitude of
% a forced linear system.
% D is 2nx2n the stiffness/damping matrix
% M is the 2nx2n mass matrix
% f is the 2n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector ‘amp’, giving the amplitude of
% each degree of freedom, and a second vector ‘phase’,
% which gives the phase of each degree of freedom

Y0 = (D+M.*1i*omega)\f;  % The i here is sqrt(-1)
% We dont need to calculate Y0bar - we can just change the sign of
% the imaginary part of Y0 using the 'conj' command
for j =1:length(f)/2
    amp(j) = sqrt(Y0(j)*conj(Y0(j)));
    phase(j) = log(conj(Y0(j))/Y0(j))/(2*i);
end

end