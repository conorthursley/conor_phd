%% generate a signal to use the bandpass filter for
% filter then makes the signal good enough to put through the APDL
% transient analysis

%highest frequency of interest we want to know
hf=50;
dt=1/(10*hf);
%length of time we want to find response for
lt=10;
steps=lt/dt;
x=linspace(0,10,steps);
%generate random numbers between [-1,1] with a length of steps
r=-1+2.*rand(steps,1);
% put it through the filter we created in bandpass filter toolbox thing
y=doFilter2(r);
% create vector to export to a csv file for APDL to read
Y=[x',y];

plot(x,r,'b',x,y,'g')

file='C:\ANSYS\Temp\Validation\DuffingValDec17\bandpassFilter.csv';
csvwrite(file,Y);


