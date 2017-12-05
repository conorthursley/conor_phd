%% solving the duffing equation theoretically to get our amplitude vs frequency curves

clear all
% close all
m=1;
kl=1e+3; %linear stiffness, k2
kn=-1; %nonlinear stiffness
c=1; %damping
f=200; %forcing amplitude 
m1=0.1;
m2=0.9;
% w=150;
mst=m1+m2;
%-------------------------------------------------
w0=sqrt(kl/m2);
A=1;
% omega=w/w0;
% syms u

% n=1;
% % for w=4:0.1:8
% w=5;
%     equ=w0+3/8*(kn/w0)*u^2+sqrt((-1+A/(u^2)))-w;  
%     dis=solve(equ,u); %,'MaxDegree',3);
%     n=n+1;
%     S=double(dis);
% % end

%-------------------------------------------------
%% Displacement parameter 
% Y=cuberoot(sqrt(X+H)+G)
row=0;
% for omega=0:0.1:3.5
%     row=row+1;
%     X=(64/(243*kn^3))*(((kl^3)/3)-((m2^3*w0^6*omega^6)/3)-kl^2*m2*w0^2*omega^2+kl*m2^2*w0^4*omega^4);
%     H=(4*m2^2*w0^4*omega^4*A^2)/(9*kn^2);
%     G=(2*m2*w0^2*omega^2*A)/(3*kn);
%     Y=G+sqrt(X+H);
%     if isreal(Y)==1;
%         V=nthroot(Y,3);
%         array(row,1)=V;
%     else 
%         % complex number roots
%         theta=atan2(imag(Y),real(Y));
%         r=sqrt(abs(Y));
%         for m=1:3
%             cis=2*pi*m/3+theta/3;
%             root1=r^(1/3)*(cos(cis)+1i*sin(cis));
%             array(row,m)=root1;
%         end
%     end
%     
%    
% end
%% Effective Mass calculations
% m_eff =[m_eff1 m_eff2]
% upsilon is the displacement parameter calculated above
%------------------------------------------------------------
% linear case

w=linspace(0,200,10000);
var3=(w/w0).^2;
omega=(w/w0);
var4=1-var3;
m_eff_linear=1+((m2/m1)/(1+(m2/m1))*((var3))./(var4));
%------------------------------------------------------------
% displacement calcs
X=(64/(243*(kn^3)))*(((kl^3)/3)-((m2^3*w0^6*omega.^6)/3)-(kl^2*m2*w0^2*var3)+(kl*m2^2*w0^4*omega.^4));
H=(4*m2^2*w0^4*omega.^4*A^2)/(9*(kn^2));
G=(2*m2*w0^2*var3*A)/(3*kn);
row=0;
for duu=1:length(w)
    row=row+1;
    S(duu)=G(duu)+sqrt(X(duu)+H(duu));
    if isreal(S(duu))==1
        V=nthroot(S(duu),3);
        array(row,1)=V;
    else
        % complex number roots
        theta=atan2(imag(S(duu)),real(S(duu)));
        r=sqrt(abs(S(duu)));
        for m=1:3
            cis=2*pi*(m/3)+(theta/3);
            root1=r^(1/3)*(cos(cis)+1i*sin(cis));
            array(row,m)=root1;
        end
    end
        
end

%------------------------------------------------------------
% new effective mass calcs
Y=array(:,1);
upsilon=Y.^2;
%------------------------------------------------------------
% First solution (no plus/minus)
m_eff1=mst-(4*m2^2*w0^2*(var4)-(9*kn*m2*upsilon'))./(4*kn*Y'*A);
% alpha1=real(m_eff1);
alpha1=imag(m_eff1);
% alpha1=abs(m_eff1);
% alpha1=[mmf1A_1;mmf1A_2;mmf1A_3];
%------------------------------------------------------------
% Second solution (no plus/minus)
Y=array(:,2);
upsilon=Y.^2;
m_eff2=mst-(4*m2^2*w0^2*(var4)-(9*kn*m2*upsilon'))./(4*kn*Y'*A);
alpha2=imag(m_eff2);
% alpha2=real(m_eff2);
% alpha2=abs(m_eff2);
% alpha2=[mmf21;mmf22;mmf23];
%------------------------------------------------------------
% Third solution (no plus/minus)
Y=array(:,3);
upsilon=Y.^2;
m_eff3=mst-(4*m2^2*w0^2*(var4)-(9*kn*m2*upsilon'))./(4*kn*Y'*A);
alpha3=imag(m_eff3);
% alpha3=real(m_eff3);
% alpha3=abs(m_eff3);
% alpha3=[mmf31;mmf32;mmf33];
%------------------------------------------------------------
% Second solution (plus/minus)
% first solution of displacement
Y=array(:,1);
upsilon=Y.^2;
m_eff2_A=mst+(4*m2^2*w0^2*(var4)*(sqrt(3i)+1)+9*kn*m2*upsilon'*(sqrt(3i)-1))./(18*kn*Y'*A);
alpha4=imag(m_eff2_A);
% alpha4=real(m_eff2_A);
% alpha4=abs(m_eff2_A);
% alpha4=[mmf2A_1; mmf2A_2; mmf2A_3];
%------------------------------------------------------------
% Second solution (plus/minus)
% second solution of displacement
Y=array(:,2);
upsilon=Y.^2;
m_eff2_B=mst+(4*m2^2*w0^2*(var4)*(sqrt(3i)-1)-9*kn*m2*upsilon'*(sqrt(3i)-1))./(18*kn*Y'*A);
alpha5=imag(m_eff2_B);
% alpha5=real(m_eff2_B);
% alpha5=abs(m_eff2_B);
% alpha5=[mmf2B_1;mmf2B_2;mmf2B_3];
%------------------------------------------------------------
Y=array(:,3);
upsilon=Y.^2;
m_eff2_C=mst+(4*m2^2*w0^2*(var4)*(sqrt(3i)+1)+9*kn*m2*upsilon'*(sqrt(3i)-1))./(18*kn*Y'*A);
alpha6=imag(m_eff2_C);
% alpha6=real(m_eff2_C);
% alpha6=abs(m_eff2_C);
% alpha6=[mmf2C_1;mmf2C_2;mmf2C_3];
%------------------------------------------------------------
% Second solution (plus/minus)
% second solution of displacement
% figure
% plot((omega),m_eff_linear,'k',omega,mmf1,omega,mmf2,omega,mmf3,'b',omega,mmf2A,omega,mmf2B,omega,mmf2C)
plot(omega,m_eff_linear,'k',omega,alpha1,omega,alpha2,omega,alpha3,omega,alpha4,omega,alpha5,omega,alpha6)
grid on
xlabel 'w/w_0'
ylabel 'm_{eff}/m_{st}'
% line([0 5],[0 0],'LineWidth',1,'Color','red')
str=sprintf('Dimensionless effective mass as a function of w/w_0 with kn=%d and A=%d',kn,A);
title(str,'FontSize',14)
axis([0 3.5 -100 100])
legend 'Effective Mass' 'cubic M_{eff} solution 1B imag' 'cubic M_{eff} solution 1B real' 'cubic M_{eff} solution 1B abs' 'cubic M_{eff} solution 1A' 'cubic M_{eff} solution 2A' 'cubic M_{eff} solution 3A' 
% 
hold on
% omega=[0:0.1:3.5];
% for ii=1:11
%     upsilon=array(ii,1);
%     
%     m_eff1=mst-(4*m2^2*w0^2*(1-omega(ii))-9*kn*m2*upsilon^2)/(4*kn*upsilon*A);
%     m_eff2=mst+(4*m2^2*w0^2*(1-omega(ii))*(sqrt(3i)+1)-9*kn*m2*upsilon^2*(sqrt(3i)+1))/(18*kn*upsilon*A);
%     m_eff(ii,:)=[m_eff1 m_eff2];
%     
% end
% hold on
% plot(omega(1:11),-m_eff(1:11,1))