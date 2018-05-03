%% Plots the transient response of a mass spring damper system excited by a harmonic signal
clear all
tic
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t_end=10000;   % t limit
t=0:dt:t_end;      % [s] time scale
t_find=6000; % the time to safely assume SS has been reached 600 seconds after initial transient begins
q=find(t==6000); p=find(t==t_end); 

omega=25;
%work period/cycle
%from 8000 seconds
t_b=8000;
t_cycle=2*pi*sqrt(1/((omega*(2*pi)))^2);
t_a=((t_b+t_cycle));
[m_min,i_min]=min(abs(t(:)-t_a));
p=i_min; q=find(t==t_b); 

%% Initial conditions: x(0) = 0, x'(0)=0
initial_x    = 0;
initial_dxdt = 0;
K=1000;
mass1=0.1;
force=1;
damp=0.0002;
z=[initial_x initial_dxdt];
w_nat=sqrt(K/mass1);

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@(t,z) rhs(t,z,omega),t,z,options);
t_new=t(q:p);
toc
%% Import comparison 
% M='U:\_PhD\Datathief\Figure3-InmanEngvib\figure3.csv';
% data=csvread(M,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'b*'); %,data(:,1),data(:,2),'r');
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
%% Calculate the PSD of the time series
FFTsize=2048;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x(:,1),hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');
%% Amplitude and resonance 
% Xamp=zeros(length(omega),1)';
% Wn=zeros(length(omega),1)';
% for jj=1:length(omega)
%     Wn(jj)=omega(jj); %omega(jj)];
%     [amp]=forced_vibration(K,m,f,Wn(jj));
%     Xamp(jj)=amp;
% end
% figure
% plot(omega/(2*pi),abs(Xamp))
% xlabel('frequency \omega'); ylabel('Amplitude x');
% title('Amplitude v Frequency ')
% % legend 'ODE45' 'DataThief'
% grid on
% set(gca,'fontsize',20) 
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
%% Work and Energy functions
%-----------Results----------------------
% extract displacement amplitudes from vector, "x"
% x = [displacement1 velo1 disp2 velo2]
% displacement
x_new=x(q:p,:);
t_new=t(q:p);
m1_disp=x_new(:,1);
% velocity 
m1_velo=x_new(:,2);
%freq
w=omega*2*pi;
%------------Work-----------------------------
% Looking at the dynamic model, apply the energy approach for KE and PE and
% find the total KE and PE of the system
% KE=KE1 + KE2
% PE=PEu1+PEu1+PE(u2-u1) (PEu1 is done twice as we have two k1 springs on
% either side of the model)
% -------------KE------------
% Kinetic Energy = 0.5*m_i*v_i^2
KE=0.5*mass1*((m1_velo).^2);
% -------------PE------------
% Potential Energy = 0.5*k_i*u_i^2
PE=0.5*K*((m1_disp).^2);

% Damping loss = pi*damping*driving_freq*amp^2
dampLoss=pi*damp*(omega*2*pi)*(max(m1_disp)).^2;

% total mechanical energy
ME=PE+KE-dampLoss;
workIn=(force*sin(omega*2*pi*t_new)).*((m1_disp));
% RMS values
y1=rms(ME);
y2=rms(workIn);

%-------Carl's method------%
% harmonic solution assumed for F and x_disp
% Input Force
In=abs((force)*sin(w*t_new).*(m1_disp));
% Output motion 
% Out=force.*max(m1_disp)*(sin(w*t_new)).^2;
Out=0.5*mass1*(m1_velo.^2) + 0.5*K*(m1_disp.^2);
% Out=0.25*(max(m1_disp))^2*(mass1*w^2+K+cos(2*w*t_new)*(mass1*w^2-K));
figure
plot(t_new,In,'b',t_new,Out,'r')
%% Work and Energy Figure
figure
plot1=plot(t_new,KE,'r',t_new,PE,'b',t_new,ME,'k',[t_b,t_a],[y1, y1],'m:',...
    t_new,workIn,'g',[t_b,t_a],[y2, y2],'c:');
xlabel('t'); ylabel('Work (J)');
set(plot1,'LineWidth',2.5)
title(['Work/Energy Calculations at ', num2str(omega), ' Hz'])
grid on
axis([t_b t_a -Inf Inf]);
legend 'KE' 'PE' 'Work Output' 'Work Output RMS' 'Work Input'  'Work Input RMS'
set(gca,'fontsize',24) 
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
    function dxdt=rhs(t,x,omega)
        mass1=0.1;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.0002;     % [Ns/m] keep as a small number to fix solver errors
        f=1;           % [N] amplitude of driving force
        w=omega;

        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1)*sin((w*2*pi)*t);

        dxdt=[dxdt_1; dxdt_2];
    end