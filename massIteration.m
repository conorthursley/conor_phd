%% Function that changes the mass ratio from 0.1 to 10. Masses that are involved are the primary and secondary mass
% of the metamaterial system, m1 and m2. Can also be used to change the spring ratios
%--------------------
% Purpose of this code is to find the best mass and spring ratios for the metamaterial system with linear springs

function massIteration()
tic
control =2; %ratio for mass ratio
% spring = 0.1:0.1:1;
y=[0;0;0;0];
ti=0; tf=1; tstep=10000; %initial, final, and stepping time
tspan = linspace(ti,tf,tstep);
figure
hold on
for ii = 1:length(control)
    mr=control(ii);
    
    [t,X]=ode45(@EOM,tspan,y);
    allx(:,:,ii)=[t X];
    u = X(:,1);
    w = X(:,3);
    
    zeta=[u w];
    [alpha, beta]=meshgrid(t,mr);
%     plot3(alpha,beta,u,  'b')
%     plot3(alpha,beta,w, 'r')
    plot(t,u,'r',t,w,'b')
end
    function dX = EOM(t, y)
        % coefficients
        k1=2e6;
        m1=1;
        k2=mr*k1;
        m2=1*m1;
        
        % define X matrix, Xdot
        x1=y(1);
        x2=y(2);
        y1=y(3);
        y2=y(4);
        %natural frequency of the system
        w0=sqrt(k2/m2);
        f=5*(2*pi);  %signal input
        
        
        % Define A matrix
        
        A=[x2; -((k1/m1)+(k2/m1))*(x1)+((k2/m1)*y1); y2; ((k2/m2)*x1)-((k2/m2)*y1)];
        
        %% Input excitation force
        
        %sawtooth with random noise
        s=sawtooth(t);
        noise=awgn(s,0.5);
        % sinusodial harmonic
        H=1*(sin(f*t));
        
        B = [0 0 H/m1 0];
        
        %% output result as the equation
        dX=A+B';
        % dy=A*y;
    end
%% graph
save test.mat allx t X alpha beta
 xlabel('Time, s')
 ylabel('Displacement, m')
%  zlabel('Spring ratio')
title(['Spring Ratio of k_r = ' num2str(control) ])
axis([0 1 -4e-4 4e-4])
 hold off
 legend 'u, mass1' 'w, mass 2'
 
grid
toc


figure
temp=abs(w)./abs(u);
plot(t,temp/1000)
end


