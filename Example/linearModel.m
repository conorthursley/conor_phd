%% Simple 1D mass-in-mass system with no damping (linear springs)
% 28/06/2017 - Conor MacDonald
clear all
close all
% Time specification
t = 0:0.01:10;

% Initial conditions
y=[0;0;0;0];



% System simulation
% [t, y] = ode45(@sys, t, y1);
[t, result] = ode45(@sys, t, y);
plot(t,y(1),t,y(2),t,y(3),t,y(4))


% Plots trajectory of displacement of 2 masses
figure
plot(t, result(:,1),t, result(:,3),t, result(:,2),t, result(:,4))
% x1=velocity of x1, x2=velocity of x2, x3=acceleration of x1,
% x4=acceleration of x2.

% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('velocity and acceleration'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('x^1_1','x^1_2','x^2_1','x^2_2')
grid on

% Plots trajectory of velocity of 2 masses

figure
plot(t, result(:,2),t, result(:,4))
%format plot
xlabel('time'); % Insert the x-axis label
ylabel('accerlation'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('x''1','x''2')
grid on

X=result(:,1);
Y = fft(X);
L=1500;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs=1000;
f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Frequency Response Function plot
% FRF = fft(output)./fft(input);
% input from H = sin(2*pi*1*t);
% output from result(:,1);
% plot(t,FRF);

%% Impulse Reaction Plot
% displacement over time [0-0.1]s
function [dA] = sys(t,y)
m1=0.1;
m2=0.05;
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
dy1=x2;
dy2=x1/m1-y1/m1;
dy3=y2;
dy4=x1/m2-y1/m1;
dy=[dy1;dy2;dy3;dy4];
B=[0 0.1*sin(5*2*pi*t) 0 0];
dA=dy+B';
end
