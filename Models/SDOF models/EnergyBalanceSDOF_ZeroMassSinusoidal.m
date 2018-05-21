%% Energy Balance 0 mass sinusoidal forcing situation 

tspan=0:0.001:pi/2; %span of t values (time)
x0=[0]; %initial displacement at t=0;

[t,X]=ode45(@(t,X) rhs(t,X), tspan, x0);  

plot(t,X,'-o')
xlabel 'time,t'
ylabel 'displacement, X'

% create for loop to iterate over the entire graph to find the area under
Fx=1000*X;
h=X;
delta=h(2)-h(1); %step value
% the curve
figure 
plot(X,Fx,'-o')
xlabel 'displacement, X'
ylabel 'Force, F'


trapint=0;  % counter
% Trap rule for loop 
for jj=1:length(h)-1
    trapint=trapint+delta*((Fx(jj)+Fx(jj+1))/2);
end
trapint  %display the number

Potential=0.5*1000.*(X).^2; %compare with Potential energy term to show 
% PE and Work are the same
Fkx=Fx*delta; %(t(2)-t(1));
figure
plot(t,Potential,t,Fkx)
legend 'potential' 'force'
%% Mass-Spring system
% The equations for the spring system
    function dxdt=rhs(t,x0)
        stiff1=1000;    % [N/m]
%         dxdt_1 = 0.2*cos(0.2*t)/stiff1;  
        dxdt_1 = 5*cos(5*t)/(stiff1);
        
        dxdt=[dxdt_1];  %return variable 
    end