%% Duffing Oscillator solved using ode45
% use ode45 to numerically solve the duffing oscillator ODE equation
% x2 + delta*x1 + beta*x0 + alpha*x0^3 = gamma*cos(omega*t)
% x2= double derivative of x wrt to t.
%----------------------------------------------

%% Coefficients 
alpha=5;
beta=1;
delta=0;
gamma=8;
omega=0.65;
tspan=[0 100];
y0=[0 0]; %initial conditions
%% ODE solver
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

tic
[t, result]=ode45(@(t,y) duffing(t, y, alpha, beta, gamma, delta, omega), tspan, y0,opts);
toc


%% Plots
figure;
plot(t,result);
xlabel('Time');
ylabel('State');
grid on
a=num2str(alpha);
str ='Forced duffing oscillator of x^{''}+\deltax1 + \alphax0 + \betax0^3 = \gammacos(wt) with \alpha=a, \beta=%d, \gamma=%d, \delta=%d, \omega=%d',alpha,beta,gamma,delta,omega;
title(str,'FontSize',12)
legend({'disp','velo'},'FontSize',14)


figure;
plot(result(:,1),result(:,2));
title(str,'FontSize',14)
xlabel('Displacement');
ylabel('Velocity');
grid on
%% *******************Acceleration**********************
%-------------u1---------------
time=t;
u1=(result(:,1))';
velocity=(result(:,2))';
nn = length(time); %Assume velocity vector is same length
ta = [time(3),time',time(nn-2)];
ve = [velocity(3),velocity,velocity(nn-2)];
t1 = ta(1:nn); t2 = ta(2:nn+1); t3 = ta(3:nn+2);
v1 = ve(1:nn); v2 = ve(2:nn+1); v3 = ve(3:nn+2);
t21 = t2-t1; t32 = t3-t2; t31 = t3-t1;
v21 = v2-v1; v32 = v3-v2;
ac_u1 = (v21./t21.*t32+v32./t32.*t21)./t31; % Approx. acceleration values
% ac_u1 = M((1:length(M)),6);
%*******************Poincare Section**********************
nnnn=size(time);
%set the index of poincare points to 1
%-----------------u1------------
np_u1=1;
for i=1:nnnn(1)
        % detect the cros-section of the trajectory with the plane y1-y2
        if (ac_u1(i)>=(2*pi)/gamma); %np_u1); 
%             (time(i)>=((2*pi)/28.69)*np_u1)
            % store detected cross-section point y1,y2 to ps1,ps2
        ps_u1(np_u1,1)=u1(i);
        ps_u1(np_u1,2)=velocity(i);
        
        % increase the index of poincare point
        np_u1=np_u1+1;
        end
end
%% plot poincare section 
figure
% plot(u1,v1,'b--')
hold on
for i=1:np_u1-1
    plot(ps_u1(i,1),ps_u1(i,2),'r.')
    % use pause to folow the plot of the poincare section
%     pause(0.05);

end
grid on
title('Poincare Section','FontSize',20)
xlabel('displacement of mass1','FontSize',20)
ylabel('velocity of mass1','FontSize',20)

%% Function

function dy = duffing(t, y, alpha, beta, gamma, delta, omega)
A=y(2);
B=gamma*cos(omega*t)-delta*y(2)-beta*y(1)-alpha*y(1).^3;
dy=[A; B];
end
