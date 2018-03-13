%% Bifurcation Code
% Taken from U:\_PhD\Books\Bifurcation Diagrams.pdf

%% Incomplete 
close all
clear all
% setup
solutionInterval = 1:1000;
rRange = 3:0.00095:4;
x=1*ones(1,length(rRange));  
index = 1;

% recurrence loop
for r = rRange 
    for n = solutionInterval
        % recurrence equation 
        x(n+1,index) = r*(x(n,index))*((1-x(n,index)));
    end
    index = index + 1;
end
% truncating the solution vector so that we only see the end behavior
cutOff = 600;
x = x(cutOff:end,:);

% %% Import comparison 
% M='U:\_PhD\Datathief\Bifurcation Map\Bifurcation_ThompsonStewart\figure9.5\figure9_5_part4.csv';
% data=csvread(M,1,0);%
%% Plot
figure
plot(rRange ,x,'k.','Marker','.','MarkerSize' ,0.1)
title('Bifurcation Diagram for ${x}_{n+1}=\mu{x}_{n}(1-{x}_{n})$','FontAngle','italic','Interpreter','Latex')
xlabel('r'); ylabel('x');
grid on
legend 'MATLAB' 'Datathief'