%% Analytical Code for AMM model

n=1; % number of metamaterial units in the chain

%% Parameters for AMM unit cell
% parameters should be passed through ODE45 to avoid reptitive changing of
% inputs
mass1=0.1;		% [kg]
mass2=mass1*0.5;
stiff1=1000;    % [N/m]
stiff2=1.5*stiff1;
damp1=0.02;     % [Ns/m] keep as a small number to fix solver errors
damp2=1*damp1;
w=10*2*pi; % driving frequency

omega_axis=logspace(-1,3,4000);

%% Ratio parameters
theta=mass2/mass1;  % mass ratio 
%---------------------------------------
w2=sqrt(stiff2/mass2);
w1=sqrt(stiff1/mass1);
%---------------------------------------
% exctiation frequency 
w=linspace(0,500,20000);
%---------------------------------------
meff=mass1+((mass2*w2^2)./(w2^2-w.^2)); % M_effective - effective mass ratio 
eta_r=w/w2;       %ratio of excitation frequency to mass2 frequency 
% this value will be the independent varaible for the bandgap diagram
eta_s=w2/w1;        % structural frequency ratio
%---------------------------------------
Qa=acos(1-(((1-eta_r.^2+theta)/(2*(1-eta_r.^2)))*eta_r.^2*eta_s^2));
qa=2*asin(sqrt((meff.*w)/(4*stiff1)));
qL=acos(1-(meff.*w.^2/(2*stiff1)));
qL(imag(qL)~=0) = nan;

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
omega_real=omega_axis/(2*pi);
plot1=loglog(omega_real,X(2*n,:),omega_real,X(1,:));
title('Amplitude Response')
grid on
xlabel('Frequency,  (Hz)')
ylabel('Amp')
set(plot1,'LineWidth',1.5)
legend 'Outer mass' 'Inner mass'
set(gca,'fontsize',20)

%-----------------Transmission

U=X((2*n),:)./(X(1,:));
%-----------------Import comparison 
N1='U:\_PhD\Datathief\Experimental_displacement_Yao\Ratio of displacements\Ratio_displacements_yao_up.csv';
data=csvread(N1,1,0);
N='U:\_PhD\Datathief\Experimental_displacement_Yao\Ratio of displacements\Ratio_displacements_yao_down.csv';
data1=csvread(N,1,0);

figure
plot2=plot(omega_axis/(2*pi),(U)); %,data(:,1),data(:,2),'g',data1(:,1),data1(:,2),'g');
title('Ratio of displacement amplitudes $|x|/|X|$','FontAngle','italic','Interpreter','Latex')
grid on
xlabel('Frequency ratio, (Hz)')
ylabel('Transmission')
legend 'MATLAB' %'Yao Experimental'

%% Dispersion Diagram
figure
plot(qL/pi,w/(2*pi)); %,w/w2,abs(imag(qL)/(pi))); 
% plot((Qa)/pi,eta_r)
title_text=['Band structures having $\theta=',num2str(theta),'$'];
title(title_text,'FontAngle','italic','Interpreter','Latex')
grid on
xlabel('Wavenumber')
ylabel('Frequency \omega')
legend_text=['\eta=',num2str(eta_s),''];
legend(legend_text,'FontAngle','italic','Interpreter','Latex')

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