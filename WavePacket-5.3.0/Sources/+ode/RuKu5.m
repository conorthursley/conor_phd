function xout = RuKu5(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Runge-Kutta 5-th order integrator
%
% The RuKu5 method is a fifth-order method, meaning that the local 
% truncation error is on the order of O(h^6), while the total accumulated 
% error is order O(h^5), https://en.wikipedia.org/wiki/Runge-Kutta_methods
%
% here: using linear approximation for u(t+tau)
%--------------------------------------------------------------------------

k1 = feval(RHS, u +   0           , xin                                                 );
k2 = feval(RHS, u + 01/04*dudt*tau, xin+tau*(01/0004)*(     k1)                         );
k3 = feval(RHS, u + 03/08*dudt*tau, xin+tau*(03/0032)*(     k1+00003*k2)                );
k4 = feval(RHS, u + 12/13*dudt*tau, xin+tau*(12/2197)*(0161*k1-00600*k2+00608*k3)       );
k5 = feval(RHS, u +   1  *dudt*tau, xin+tau*(01/4104)*(8341*k1-32832*k2+29440*k3-845*k4));

xout = xin + tau * ((25/216)*k1+(1408/2565)*k3+(2197/4104)*k4-(1/5)*k5);

end
