function xout = Heun(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Midpoint explicit integrator
% https://en.wikipedia.org/wiki/Heun%27s_method
%--------------------------------------------------------------------------

k1 = feval(RHS, u,          xin       );
k2 = feval(RHS, u+dudt*tau, xin+k1*tau);

xout = xin + tau * (k1 + k2) / 2;

end
