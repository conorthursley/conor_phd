function xout = RuKu4(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Runge-Kutta 4-th order integrator
%
% The RuKu4 method is a fourth-order method, meaning that the local 
% truncation error is on the order of O(h^5), while the total accumulated 
% error is order O(h^4), https://en.wikipedia.org/wiki/Runge-Kutta_methods
%
% here: using linear approximation for u(t+tau)
%--------------------------------------------------------------------------

k1 = feval(RHS, u +  0          , xin +  0        ); % original RuKu: t
k2 = feval(RHS, u + 1/2*dudt*tau, xin + 1/2*k1*tau); % original RuKu: t+tau/2
k3 = feval(RHS, u + 1/2*dudt*tau, xin + 1/2*k2*tau); % original RuKu: t+tau/2
k4 = feval(RHS, u +  1 *dudt*tau, xin +  1 *k3*tau); % original RuKu: t+tau

xout = xin + tau * (k1 + 2*k2 + 2*k3 + k4)/6;

end
