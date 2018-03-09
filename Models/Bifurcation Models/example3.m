%% Bifurcation Code
% Taken from U:\_PhD\Books\Bifurcation Diagrams.pdf
close all
clear all
% setup
solutionInterval = 1:1000;
rRange = 2.4:0.00095:4.0;
x=0.585*ones(1,length(rRange));  
index = 1;

% recurrence loop
for r = rRange 
    for n = solutionInterval
        x(n+1,index) = r*x(n,index)*(1-x(n,index));
    end
    index = index + 1;
end
% truncating the solution vector so that we only see the end behavior
cutOff = 600;
x = x(cutOff:end,:);

%% Import comparison 
M='U:\_PhD\Datathief\Bifurcation Map\Bifur_4.csv';
data=csvread(M,1,0);
%% Plot
figure
plot(rRange ,x,'k.','Marker','.','MarkerSize' ,0.1)