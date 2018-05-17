%% Energy Balance finite mass situation 

tspan=0:0.01:1; %span of t values (time)
x0=[0]; %initial displacement at t=0;

[t,X]=ode45(@(t,X) rhs(t,X), tspan, x0);  

plot(t,X,'-o')

% create for loop to iterate over the entire graph to find the area under
% the curve

h=t;
delta=h(2)-h(1); %step value
trapint=0;  % counter
% Trap rule for loop 
for jj=1:length(h)-1
    trapint=trapint+delta*((X(jj)+X(jj+1))/2);
end
trapint  %display the number

Potential=0.5*1000.*X(end).^2 %compare with Potential energy term to show 
% PE and Work are the same

%% Mass-Spring system
% The equations for the spring system
    function dxdt=rhs(t,x0)
        stiff1=1000;    % [N/m]
        forceAmp=1;    % forcing amplitude
        
        dxdt_1 = -(forceAmp*cos(t))/stiff1;  %*** F=Asint
%         dxdt_2 = -(forceAmp*sin(t))/stiff1;
        
        dxdt=[dxdt_1];  %return variable 
    end