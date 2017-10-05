 %% Function that changes the mass ratio from 0.1 to 10. Masses that are involved are the primary and secondary mass
 % of the metamaterial system, m1 and m2. Can also be used to change the spring ratios
 %--------------------
 % Purpose of this code is to find the best mass and spring ratios for the metamaterial system with linear springs
 
 function massIteration()
 tic
 control = 0.1:0.1:1; %ratio for mass ratio
 spring = 0.1:0.1:1;
 y=[0;0;0;0];
 ti=0; tf=10; tstep=100; %initial, final, and stepping time
 tspan = linspace(ti,tf,tstep);
 hold on
 for ii = 1:length(control) 
  
         
         [t,X]=ode45(@(t,X)EOM(t,y,control(ii)),tspan,y);
         allx(:,:,ii)=[t X];
         u = X(:,1);
        
         
     
     w = X(:,3);
     plot3(t,u,w)
 end
     function dX = EOM(t, y, mr)
    % coefficients
        k1=2e6;
        m1=1;
        k2=0.1*k1;
        m2=mr*m1;
    
    % define X matrix, Xdot
        x1=y(1);
        x2=y(2);
        y1=y(3);
        y2=y(4);
    %natural frequency of the system
        w0=sqrt(k2/m2);
        f=5/(2*pi);  %signal input
       
      
    % Define A matrix
        
        A=[x2; -((k1/m1)+(k2/m1))*(x1)+((k2/m1)*y1); y2; ((k2/m2)*x1)-((k2/m2)*y1)];
    
    %% Input excitation force
    
    %sawtooth with random noise
        s=sawtooth(t);
        noise=awgn(s,0.5);
% sinusodial harmonic
        H=2*(sin(f*t));
    
        B = [0 0 H/m1 0];
    
    %% output result as the equation
        dX=A+B';
    % dy=A*y;
     end

 save test.mat allx w u t X
%  xlabel('Time, s')
%  ylabel('Displacement, m')
%  hold off
%  legend 'u, mass1' 'w, mass 2'
 grid
 toc
 end

 
