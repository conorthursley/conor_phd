%% generate a signal to use the bandpass filter for
% filter then makes the signal good enough to put through the APDL
% transient analysis
clear all
%highest frequency of interest we want to know
hf=100;
dt=1/(1000*hf);
%length of time we want to find response for
lt=5;
steps=(lt/dt);
x=linspace(0,lt,(steps+1));
% x=0:1/steps:lt;
%generate random numbers between [-1,1] with a length of steps
A=20;
r=-A+(2*A).*rand(5e5,1);
% put it through the filter we created in bandpass filter toolbox thing
% y=doFilter2(r,dt);
% create vector to export to a csv file for APDL to read
Y=[x',r];

plot(x,r,'b') %,x,y,'g')

file='C:\ANSYS\Temp\Validation\DuffingValDec17\HigherAmp.csv';
csvwrite(file,Y);


