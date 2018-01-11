%% FEM code for metamaterial units
%---------------------------------------------------
% n = no.of unit cells.
% each unit cell contains 2 masses
n=2; % must be greater than 2
%---------------------------------------------------
% Type of metamaterial system
dyn=0; % read before the ODE solver strarts to determine which stiffness matrix to use in the function
% 0 is a linear system, 1 is a piecewise system, 2 is a nonlinear system
%---------------------------------------------------
% Time span
tic
tspan=[0 10];
%---------------------------------------------------
%% set parameters for stiffness, mass and IC
m1=0.1;
m2=0.5*m1;  %kg
k1=1000;
k2=320;  % N/m
y0=zeros(1,(n*2));
% y0(3)=1e-3; %initial displacement of first mass
%---------------------------------------------------
% Natural frequencies for inner and outer masses
w_1=sqrt(k1/m1);
w_2=sqrt(k2/m2);
%---------------------------------------------------
% Input frequency, Hz and amplitude
f_ratio=1;
freq = (w_2/2*pi)*f_ratio; %Hz
A=0.1; % Amplitude
input = [freq A];
%---------------------------------------------------
% Mass matrix is computed in a function at the bottom of the code
%---------------------------------------------------
%% Stiffness matrix
% creat with n*2 (n is no. of unit cells, with 2 masses and 2 springs in
% each cell)
K=zeros(n*2);
% function populates the stiffness matrix based on the spring settings in
% the ODE function (nonlinear, piecewise, linear, etc)
%---------------------------------------------------
%% ODE options
opts=odeset('Mass',m(m1,m2,n),'OutputFcn',@odeplot); %, 'Events', @events);
%---------------------------------------------------
%% ODE initialisation 
% switch cases for liner, piecewise, and nonlinear systems
switch dyn
    case 0
        disp('Linear case')
        %---------------------------------------------------
        % Determine Stiffness Matrix before calling ODE function
        % every first row follows: k2(u1-u2) + k1(2*u1-u1(cell to the left)
        % -u1(cell to the right))
        % every second row follows: k2(u2-u1)
        % except the first row and last 2x2 diagonal, iteration will fill the
        % matrix
        %---------------------------------------------------
        % first row, uses the first cell's secondary mass's displacement, y(2)
        %calculate the stiffness using the piecewise function
        % populate
        K(1,1)=k2+(k1);
        K(1,2)=-k2;
        K(1,3)=-k1;
        %---------------------------------------------------
        % last 2x2 diagonal, uses the last cell's secondary mass, y(end)
        % populate
        q=2*n;
        K(q,q)=k2;
        K(q-1,q-1)=k2+k1;
        K(q-1,q)=-k2;
        K(q,q-1)=-k2;
        %---------------------------------------------------
        % Iteration of remainder of cell, ignoring the first row and last diagonal
        for i=2:((2*n)-2)
            if mod(i,2)==1  % odd number and therefore a primary mass
                K(i,i)=k2+(2*k1);
                K(i,i+1)=-k2;
                K(i,i+2)=-k1;
            elseif mod(i,2)==0 % even number and therefore a secondary mass
                K(i,i)=k2;
                K(i,i-1)=-k2;
                K(i+1,i-1)=-k1;
            end
        end
        %---------------------------------------------------
        % Call the ODE solver
        
        for jj=1:50
            freq=jj/10;
            input=[freq A];
            
            [t,y] = ode45(@(t,y) L_f(t,y,K,n,input),tspan,y0,opts); %te,ye,ie
            x=y(:,(end-1));
            
            pks=findpeaks(x);
            xmax_norm(jj,:)=mean(pks)/A;
            omega_norm(jj,:)=(freq*2*pi)/w_2;
            
        end
        figure
        plot(omega_norm,xmax_norm,'ro-')
        grid on
        %---------------------------------------------------
        toc
    case 1
        disp('Piecewise case')
        [t,y,te,ye,ie] = ode45(@(t,y) P_f(t,y,K,k1,k2,n,input),tspan,y0,opts);
        toc
    case 2
        disp('Non-linear case')
        [t,y,te,ye,ie] = ode45(@(t,y) f(t,y,K,k1,k2,n,input),tspan,y0,opts);
        toc
    otherwise 
        error('Error, check dynamic model assignment number, dyn')
end

%---------------------------------------------------
%% Figures/plots
% Displacement
figure
plot(t,y(:,1)/A,'r-',t,y(:,2)/A,'b-s',t,y(:,3)/A,'g-',t,y(:,4)/A,'k-s')
grid on
legend 'mass_1' 'mass_2' 'mass_3' 'mass_4' 'mass_5'
%---------------------------------------------------
% Frequency Response


%---------------------------------------------------
% Snapshot
% snapshot of the chain at time, t
% point=1:1:n; % discretise n into the number of cells for the plot
% figure
% plot(point,y(end,point))
% grid on
% legend 'displacement of chain' 
%---------------------------------------------------
% FFT
% Change the t in the below line to be either a loop or some workaround.
% error is getting thrown as t is a large vector. for loop might be ideal. 
% h=input(2)*(sin(2*pi*input(1)*tspan+0.2))+2*input(1)*sin(2*pi*input(1)*2*tspan+2.2)+rand(1,tspan);
% H=fft(h);
% H_mag=abs(H);
% dt=mean(diff(t));  %average time step done in the ode45 computation
% Fs=1/dt;
% n=length(t);  %length of signal = number of samples
% q=pow2(nextpow2(n));  %transform length
% dft=fft(h,q)/n; % DFT of signal
% fr = (0:q-1)*(Fs/q)/10;
% fourier = abs(dft);     
% figure
% plot(fr(1:floor(q/2)),fourier(1:floor(q/2)))
% title('Single-Sided Amplitude Spectrum of U1(t)')
% grid on
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

%% ODE integrator functions
% Case 0 - Linear Function
% Case 1 - Piecewise Function
% Case 2 - Non-linear Function
%% Case 0 - Linear ODE Function
function dy = L_f(t,y,K,n,input)

%---------------------------------------------------
% Create the signal input to act on the last cell's primary mass
B=zeros(1,(2*n));
B((2*n)-1)=input(2)*(sin(2*pi*input(1)*t)); %+0.2)+sin(2*pi*input(1)*2*t+2.2));
%---------------------------------------------------
dy = -K*y+B';
end
%% Case 1 - ODE Piecewise Function
function dy = P_f(t,y,K,k1,k2,n,input)
%---------------------------------------------------
% Stifness calculator for piecewise function for secondary mass springs, k2
syms Q
%---------------------------------------------------
% Piecewise case
spring=piecewise(Q<=5e-4,k1,5e-4<Q<1.5e-3, k2, Q>=1.5e-3, k1); 
%---------------------------------------------------
% each time the iteration needs to determine the spring value, the loop
% calls upon the above systematic funcion to do so.
%---------------------------------------------------
% every first row follows: k2(u1-u2) + k1(2*u1-u1(cell to the left)
% -u1(cell to the right))
% every second row follows: k2(u2-u1)
% except the first row and last 2x2 diagonal, iteration will fill the
% matrix
%---------------------------------------------------
% first row, uses the first cell's secondary mass's displacement, y(2)
k=subs(spring, Q, y(2));
k_one=double(k); %calculate the stiffness using the piecewise function
% populate
K(1,1)=k_one+(2*k1); 
K(1,2)=-k_one;
K(1,3)=-k1;
%---------------------------------------------------
% last 2x2 diagonal, uses the last cell's secondary mass, y(end)
k=subs(spring, Q, y(end));
k_n=double(k);
% populate
q=2*n;
K(q,q)=k_n;
K(q-1,q-1)=k_n+k1;
K(q-1,q)=-k_n;
K(q,q-1)=-k_n;
%---------------------------------------------------
% Iteration of remainder of cell, ignoring the first row and last diagonal
for i=2:((2*n)-2)
    if mod(i,2)==1  % odd number and therefore a primary mass
        j=i+1; % the next secondary mass that is coupled to the primary
        k=subs(spring, Q, y(j+1));
        k2=double(k);
        K(i,i)=k2+(2*k1);
        K(i,i+1)=-k2;
        K(i,i+2)=-k1;
    elseif mod(i,2)==0 % even number and therefore a secondary mass
        k=subs(spring, Q, y(i));
        k2=double(k);
        K(i,i)=k2;
        K(i,i-1)=-k2;
        K(i+1,i-1)=-k1;
    end
end
%---------------------------------------------------
% Create the signal input to act on the last cell's primary mass
B=zeros(1,(n*2));
B((2*n)-1)=heaviside(input(2)*(sin(2*pi*input(1)*t)+sin(2*pi*6*input(1)*t)+sin(2*pi*4*input(1)*t)));
%---------------------------------------------------
dy = -K*y+B';
end

%% Mass function

function M = m(m1,m2,n)
%---------------------------------------------------
% Mass matrix
m_vector=[m1 m2];
m_array=repmat(m_vector,n);
M=diag(m_array(1,:));
end

%% Event function 

% function [value,isterminal,direction] = events(~,y)
%       value = [double(abs((y(n)-(y(n-1))))<=(5e-4))];
%            % detect when the bounds gets crossed
%       isterminal = [0]; % halt integration, reverse direction 
%       % 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
%       direction = [0]; % approaching the event from any which way
%       %0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing, ...
%       % ... and -1 locates only zeros where the event function is decreasing
% end

