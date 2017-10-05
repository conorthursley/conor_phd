 %% Function that changes the frequency to find the amplitude ratio 
 %--------------------
 % Purpose of this code is to prepare for f,m, and k optimisation
 % x/X=k2/(k2-m2*w^2)
 
 function FrequencyIteration()
 tic
 f1=1; fstep=0.5; f2=50; %initial, stepping, final frequency 
 control = f1:fstep:f2; % frequency range from 0 to 200 hz (low frequency range)
 hold on
 for ii = 1:length(control) 
    y=[0;0;0;0];
    ti=0; tf=10; tstep=1000; %initial, final, and stepping time
    tspan = linspace(ti,tf,tstep);
    freq=control(ii);
    [t,X] = ode45(@EOM,tspan,y);
    freqx(ii,:,:)=[t X];

    
    u = X(:,1);
    w = X(:,3);
    fax=linspace(f1,f2,length(t)); %frequency axis, fax
    plot(fax,u,'b',fax,w,'r')

 end
     function dX = EOM(t, y)
    % coefficients
        k1=2e6;
        m1=1;
        k2=0.1*k1;
        m2=0.9*m1;
    
    % define X matrix, Xdot
        x1=y(1);
        x2=y(2);
        y1=y(3);
        y2=y(4);
    %natural frequency of the system
        w0=sqrt(k2/m2);
        f=freq/(2*pi);  %signal input
       
      
    % Define A matrix
        
        A=[x2; -((k1/m1)+(k2/m1))*(x1)+((k2/m1)*y1); y2; ((k2/m2)*x1)-((k2/m2)*y1)];
    
%% Input excitation force
    
    %sawtooth with random noise
%         s=sawtooth(t);
%         noise=awgn(s,0.5);
     %sinusodial harmonic   
        H=1*(sin(f*t));
    
        B = [0 0 H/m1 0];
    
%% output result as the equation
        dX=A+B';
    % dy=A*y;
     end

 save freq.mat freqx 
 xlabel('Frequency, Hz')
 ylabel('Displacement, m')
 hold off
 legend 'u, mass1' 'w, mass 2'
 grid
 toc
 end

 
