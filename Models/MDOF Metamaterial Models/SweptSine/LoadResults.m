%% Load Results and plot the saved data
load 5UweaklyNL.mat

%% Model Parameters
mass1=0.1;		% [kg]
mass2=0.5*mass1;
stiff1=1000;    % [N/m]
stiff2=1500;

w2=sqrt(stiff2/mass2)/(2*pi);
theta=mass2/mass1;

%% Results
% extract displacement amplitudes from vector, "amplitude"
% amplitude = [displacement1 velo1 disp2 velo2 disp'n' velo'n' .....]
%----------First unit cell--------------
% displacement
m1_disp=amplitude(:,1);
m2_disp=amplitude(:,3);
% velocity 
m1_velo=amplitude(:,2);
m2_velo=amplitude(:,4);
%----------Last unit cell---------------
% displacement, 'last mass' - 'lm'
lm1_disp=amplitude(:,end-3);
lm2_disp=amplitude(:,end-1);
% velocity 
lm1_velo=amplitude(:,end-2);
lm2_velo=amplitude(:,end);
%% Plot the results
% Plot the time series
figure
plot1=plot(swept_sine_range,mag2db(m1_disp),swept_sine_range,mag2db(m2_disp)...
    ,swept_sine_range,mag2db(lm1_disp),swept_sine_range,mag2db(lm2_disp));
set(plot1,'LineWidth',2)
xlabel('Freqency range (Hz)'); ylabel('Transmittance (dB)');
title('Swept Sine results for displacement')
grid on
legend 'mass1' 'mass2' 'last mass1' 'last mass2'
set(gca,'fontsize',20) 
%% Non-dimensionalised Results
% non-dimensionalise the frequency 
NDf=swept_sine_range/w2;
% Plot the time series
figure
plot1=plot(NDf,mag2db(m1_disp),NDf,mag2db(m2_disp)...
    ,NDf,mag2db(lm1_disp),NDf,mag2db(lm2_disp));
set(plot1,'LineWidth',2)
xlabel('\eta'); ylabel('Transmittance (dB)');
title('Transmittance curve for weakly NL MIM model with NL = 1600')
grid on
legend 'mass1' 'mass2' 'last mass1' 'last mass2'
set(gca,'fontsize',24) 

%% Transmission
% U=Xn/X1
U=abs(lm2_disp)./(m2_disp);
figure
plot(swept_sine_range,20*log10((U)),'b')
xlabel('Freqency range (Hz)'); ylabel('Transmittance (dB)');
title('Transmittance of the finite system')
grid on
set(gca,'fontsize',20) 
% periodogram(x(:,3))
%% Plot the Kinetic Energy ratio
%------------------numerical
KE1=0.5*mass1.*(m1_velo.^2);
KE2=0.5*mass2.*(m2_velo.^2);
RDR=KE2./(KE1+KE2);   %ratio
%------------------analytical
RDR_a=theta./((1-(swept_sine_range./w2).^2).^2+theta);

figure
plot1=plot(swept_sine_range/w2,RDR); hold on
plot2=plot(swept_sine_range/w2,RDR_a,'r--');
set([plot1 plot2],'LineWidth',3.5)
xlabel('Normalised freqency \omega/\omega_2'); ylabel('Ratio of Kinetic Energy');
title('Swept Sine results for energy distribution rate of the AMM model')
grid on
legend_text=['\theta=',num2str(theta)];
legend(legend_text,'Analytical Result','FontAngle','italic','Interpreter','Latex')
set(gca,'fontsize',20) 