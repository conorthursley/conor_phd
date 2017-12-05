%
% solve the DBR photonic band structure for a simple
% 1D DBR. air-spacing d, periodicity a, i.e, a > d,
% we assume an infinite stack of 1D alternating eps_r|air layers
% y-polarized, z-directed plane wave incident on the stack
% periodic in the z-direction; 
%

%parameters
d=8; %air gap
a=10; %total periodicity
d_over_a = d/a;
eps_r =12.2500; %dielectric constant, like GaAs,

% max F.S coefs for representing E field, and Eps(r), are
Mmax=50;

% Q matrix is non-symmetric in this case, Qij != Qji
% Qmn = (2*pi*n + Kz)^2*Km-n
% Kn = delta_n / eps_r + (1 - 1/eps_r)(d/a)sinc(pi.n.d/a)
% here n runs from -Mmax to + Mmax,

freqs=[];
for Kz=-pi/a:pi/(10*a):+pi/a
    Q=zeros(2*Mmax + 1);
    for x=1:2*Mmax+1
        for y=1:2*Mmax+1
            X=x-Mmax;
            Y=y-Mmax;
            kn=(1 -1/eps_r)*d_over_a.*sinc((X-Y).*d_over_a) + ((X-Y)==0)*1/eps_r;
            Q(x,y)=(2*pi*(Y-1)/a + Kz).^2*kn;% -Mmax<=(Y-1)<=Mmax
        end
    end
    
    fprintf('Kz = %g\n',Kz)
    omega_c=eig(Q);
    omega_c=sort(sqrt(omega_c));%important step.
    freqs=[freqs; omega_c.'];
end

close()
figure()
hold on
idx=1;

for idx=1:length(-pi/a:pi/(10*a):+pi/a)
    plot(-pi/a:pi/(10*a):+pi/a,freqs(:,idx),'.-')
end
    
hold off
xlabel('Kz')
ylabel('omega/c')
title(sprintf('PBG of 1D DBR with d/a=%g, Epsr=%g',d/a,eps_r))