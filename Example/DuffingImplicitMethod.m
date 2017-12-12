%% plotting duffing relation
% 
close all

%parameters
a=1;
y=3/4*0.5;
d=0;
b=[-0.5 -0.05 0 0.05 0.05]; %linear oscillator - as B is the nonlinear term
% z is steady state displacment and w is frequency
figure
for i=1:5
   B=b(i);
   
   f = @(w,z) ((w.^2-a-(3/4)*B.*z.^2)^2+(d.*w).^2).*z.^2-y^2;

   fimplicit(f,[0 3.5 0 15])
   hold on
end

xlabel('frequency, w','FontSize',14) % x-axis label
ylabel('amplitude, z','FontSize',14) % y-axis label
title('Duffing Relation $\ddot{x}$ + $\delta\dot{x}$ + $\alpha{x}$ $\beta{x}^3$ = $\gamma{cos(\omega{t})}$', 'Interpreter','latex' ,'FontSize',14)
legend({'\beta = %d','\beta = %d','\beta = %d','\beta = %d','\beta = %d'},'FontSize',14)
grid on
