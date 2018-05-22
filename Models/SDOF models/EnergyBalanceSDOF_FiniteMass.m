%% Energy Balance finite mass situation 

tspan=0:0.001:1; %span of t values (time)
x0=[0,0]; %initial displacement at t=0;

[t,X]=ode45(@(t,X) rhs(t,X), tspan, x0);  
xdot=rhs(t,X.').';  %acceleration
plot(t,X,'-o') %,t,xdot(:,2))

% create for loop to iterate over the entire graph to find the area under
% the curve
Potential=0.5*1000.*(X(:,1)).^2; %compare with Potential energy term to show 
% PE and Work are the same
Kinetic=0.5*0.1.*(X(:,2)).^2;
Fx=Potential+Kinetic;
plot(t,Fx,t,Potential,t,Kinetic)

h=t;
delta=h(2)-h(1); %step value
trapint=0;  % counter
% Trap rule for loop 
for jj=1:length(h)-1
    trapint=trapint+delta*((Fx(jj)+Fx(jj+1))/2);
end
trapint  %display the number


%% Mass-Spring system
% The equations for the spring system
    function dxdt=rhs(t,x0)
        stiff1=1000;    % [N/m]
        forceAmp=1;    % forcing amplitude
%         damp1=0;
        mass1=0.1;
%         w=0.5;
        
        dxdt_1 = x0(2,:);  %*** F=Asint
        dxdt_2 = - (stiff1/mass1)*x0(1,:) +...
            (forceAmp/mass1);
        
        dxdt=[dxdt_1;dxdt_2];  %return variable 
    end