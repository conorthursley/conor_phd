%% generate a signal to use the bandpass filter for
% filter then makes the signal good enough to put through the APDL
% transient analysis
clear all
%highest frequency of interest we want to know
hf=100;
dt=1/(10*hf);
%length of time we want to find response for
lt=500;
steps=(lt/dt);
x=linspace(0,lt,(steps));
% x=0:1/steps:lt;
%generate random numbers between [-1,1] with a length of steps
A=1;
% r=-A+(2*A).*rand(steps,1);
r=A.*sin(16.3611*2*pi*x);
% r = Amplitude*sin(w*t), where t is the timesteps, or X

% put it through the filter we created in bandpass filter toolbox thing
% y=doFilter2(r,dt);
% create vector to export to a csv file for APDL to read
Y=[x' r'];
figure
plot(x,r,'b') %,x,y,'g')

file='U:\_PhD\APDL\Validation\DuffingValDec17\NLmodeC4.csv';
csvwrite(file,Y);


