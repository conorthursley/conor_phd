function xout = Euler1(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Euler 1-st order explicit integrator
% https://en.wikipedia.org/wiki/Euler_method
%--------------------------------------------------------------------------

k = feval(RHS, u, xin);
xout = xin + tau * k;

end
