%% Plots the transient response of a mass spring damper system excited by a harmonic signal
clear all
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:1000;      % [s] time scale
omega=logspace(-1,3,length(t));
%% Initial conditions: x(0) = 0, x'(0)=0
initial_x    = 0;
initial_dxdt = 0.2;
K=1000;
m=10;
f=23;
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt],options );
%% Import comparison 
M='U:\_PhD\Datathief\Figure3-InmanEngvib\figure3.csv';
data=csvread(M,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'b*',data(:,1),data(:,2),'r');
xlabel('t'); ylabel('x');
set(plot1,'LineWidth',2)
title('Time Series')
legend 'ODE45' 'DataThief'
grid on
set(gca,'fontsize',20) 
% Plot the phase plane
figure
plot2=plot(x(:,1),x(:,2)); %,data(:,1),data(:,2),x1(:,1),x1(:,2));
set(plot2,'LineWidth',2)
xlabel('$x$','FontAngle','italic','Interpreter','Latex'); ylabel('$\dot{x}$','FontAngle','italic','Interpreter','Latex');
title('Phase potrait')
grid on
legend 'ODE45' 
set(gca,'fontsize',20) 
%%% Calculate the PSD of the time series
FFTsize=2048;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x(:,1),hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');
%% Amplitude and resonance 
Xamp=zeros(length(omega),1)';
Wn=zeros(length(omega),1)';
for jj=1:length(omega)
    Wn(jj)=omega(jj); %omega(jj)];
    [amp]=forced_vibration(K,m,f,Wn(jj));
    Xamp(jj)=amp;
end
figure
plot(omega/(2*pi),abs(Xamp))
xlabel('frequency \omega'); ylabel('Amplitude x');
title('Amplitude v Frequency ')
% legend 'ODE45' 'DataThief'
grid on
set(gca,'fontsize',20) 
%% FFT
figure
dt=abs(mean(diff(t)));  %average time step done 
Fs=1/(dt);
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1); 
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
% plot(flop,abs(y),'LineWidth',2)
title('FFT')
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)
%% *******************Poincare Section**********************
nnnn=length(t);
%set the index of poincare points to 1
%-----------------u1------------
np_u1=1;
w_f=20;      % [rad/s] - forcing frequency 
for i=1:nnnn
        % detect the cros-section of the trajectory with the plane y1-y2
        if (t(i)>=(1/(w_f/(2*pi)))*np_u1)
%             (ac_u1(i)>=(2*pi)*np_u1)
            % store detected cross-section point y1,y2 to ps1,ps2
        ps_u1(np_u1,1)=x(i,1);
        ps_u1(np_u1,2)=x(i,2);
        
        % increase the index of poincare point
                np_u1=np_u1+1;
        end
end
%% plot poincare section 

figure
hold on
plot(ps_u1(:,1),ps_u1(:,2),'r.')

% Looping the plotting points for 'effect'
% for i=1:np_u1-1
%     x_dot=ps_u1(i,2);
%     x_disp=ps_u1(i,1);
%     plot(x_disp,x_dot,'r+')
% %     hold on
%     % use pause to folow the plot of the poincare section
% %     pause(0.05);
% end
grid on
title('Poincare Section','FontSize',20)
xlabel('${x(t)}$','FontAngle','italic','Interpreter','Latex')
ylabel('$\frac{d x(t)}{dt}$','FontAngle','italic','Interpreter','Latex')
toc
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=10;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.00000001;     % [Ns/m] keep as a small number to fix solver errors
        f=23;            % [N] amplitude of driving force

        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1)*cos(2*sqrt(stiff1/mass1)*t);

        dxdt=[dxdt_1; dxdt_2];
    end
%% Force Vibration Amplitude function
function X = forced_vibration(K,M,f,omega)
% Function to calculate steady state amplitude of
% a forced linear system.
% K is nxn the stiffness matrix
% M is the nxn mass matrix
% f is the n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector X, giving the amplitude of
% each degree of freedom
%
X = (K-M*omega^2)\f;

end