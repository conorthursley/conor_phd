function xout = RuKu6(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Runge-Kutta 6-th order integrator
%
% The RuKu6 method is a sixth-order method, meaning that the local 
% truncation error is on the order of O(h^7), while the total accumulated 
% error is order O(h^6), https://en.wikipedia.org/wiki/Runge-Kutta_methods
%
% here: using linear approximation for u(t+tau)
%--------------------------------------------------------------------------

k1 = feval(RHS, u +   0           , xin                                                               );
k2 = feval(RHS, u + 01/04*dudt*tau, xin+tau*(01/0004)*k1                                              );
k3 = feval(RHS, u + 03/08*dudt*tau, xin+tau*(03/0032)*(k1+3*k2)                                       );
k4 = feval(RHS, u + 12/13*dudt*tau, xin+tau*(12/2197)*(161*k1-600*k2+608*k3)                          );
k5 = feval(RHS, u +   1  *dudt*tau, xin+tau*(01/4104)*(8341*k1-32832*k2+29440*k3-845*k4)              );
k6 = feval(RHS, u + 01/02*dudt*tau, xin+tau*(-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5));

xout = xin + tau * 1/5*((16/27)*k1+(6656/2565)*k3+(28561/11286)*k4-(9/10)*k5+(2/11)*k6);

end
