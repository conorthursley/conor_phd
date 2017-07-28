function dy=sys(t,y)
%matrix A has the spring and mass values
k1=1; 
m1=1;
k2=2*k1;
m2=9*m1;

A=[0 0 1 0;0 0 0 1; -((k1/m1)+(k2/m1)) k2/m1 0 0; k2/m2 -k2/m2 0 0];
%input excitation force
H = sin(2*pi*50*t);
%input force applied to the first mass, x1
B=[0 0 H 0];
%output result as the equation 
dy=A*y+B';
% dy=A*y;

end



