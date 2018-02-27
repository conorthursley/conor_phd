%% Function that varys the initial conditions to generate a bifurcation diagram
%
%
%--------------------

tic

tspan = [0 1];
L=0:1e-4:3e-3; %our parameter for the initial conditions
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);


parfor ii = 1:length(L)
    y=[L(ii);0;L(ii);0];
    [t,X]=ode45(@EOM,tspan,y,opts);
%     allx(:,:,ii)=[t X];
    
    for c=1:5
        g=c/5;
        [p, q]=min(abs(t-g));
        cell1(c,1,ii) = X(q,1);
        cell2(c,1,ii) = X(q,3);
    end
    
    
end
toc


%% graph
figure
hold on
for ii=1:length(L)
    temp1 = cell1(:,1,ii);
    var1 = cell2(:,1,ii);
    plot(L(ii),temp1,'r+',L(ii),var1,'b*')
end
    
xlabel('Varying parameter - Initial displacement, x')
ylabel('Displacement, m')
%  zlabel('Spring ratio')
title('Bifurcation Diagram')
hold off
legend 'u, mass1' 'w, mass 2'
 
grid on
toc

%% Function, Equations of Motion (EOM)
function dX = EOM(t, y)
% Define parameters of spring mass model
k1=2e6;
m1=1;
k2=2e5;
m2=2*m1;

%define X matrix, Xdot
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);

syms q
spring=piecewise(q<=5e-4,k1,5e-4<q<1.5e-3, k2, q>=1.5e-3, k1);
k=subs(spring, q, y1);
k3=double(k);

% w0=sqrt(k2/m2);

f=5*(2*pi); %5Hz
% Define A matrix
A=[x2; ...
    -((k1/m1)+(k3/m1))*(x1)+((k3/m1)*y1); ...
    y2; ...
    ((k3/m2)*x1)-((k3/m2)*y1)];

% Input excitation force
%sinusodial harmonic
H=1*(sin(f*t));

B = [0 0 H/m1 0];


%output result as the equation 
dX=A+B';
% dy=A*y;
end
