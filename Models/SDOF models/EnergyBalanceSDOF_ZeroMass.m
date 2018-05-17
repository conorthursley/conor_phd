%% Energy Balance 0 mass situation 

xspan=0:0.01:1; %span of x values (displacement)
F0=0; %initial forcing at x=0;

[x,F]=ode45(@(x,F) rhs, xspan, F0);  

plot(x,F,'-o')

% create for loop to iterate over the entire graph to find the area under
% the curve

h=x;
delta=h(2)-h(1); %step value
trapint=0;  % counter
% Trap rule for loop 
for jj=1:length(h)-1
    trapint=trapint+delta*((F(jj)+F(jj+1))/2);
end
trapint  %display the number

Potential=0.5*1000.*x(end).^2 %compare with Potential energy term to show 
% PE and Work are the same

%% Mass-Spring system
% The equations for the spring system
    function dxdt=rhs(x,F0)
        stiff1=1000;    % [N/m]
        dxdt_1 = stiff1;  %*** F=k*x so F'=k ***
        
        dxdt=[dxdt_1];  %return variable 
    end