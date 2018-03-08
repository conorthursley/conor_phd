close all
clear all
solutionInterval = 1:34;
x(1) = 0.5;
r = 0.8;
for n = solutionInterval 
x(n+1) = -r*x(n) + x(n)^3;
end
plot (x,'x')