close all
clear all
solutionInterval = 1:1000;
rRange = 0.8:0.001:2.6;
x = 0.5*ones (1, length ( rRange ));
index = 1;
for r = rRange 
for n = solutionInterval 
x(n+1, index ) = -r*x(n, index ) + x(n, index )^3;
end
index = index + 1;
end
cutOff = 600;
x = x( cutOff :end ,:);
figure
plot(rRange ,x,'k.','Marker','.','MarkerSize' ,0.5)