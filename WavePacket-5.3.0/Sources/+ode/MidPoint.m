function xout = MidPoint(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Midpoint explicit integrator
% see https://en.wikipedia.org/wiki/Midpoint_method
%--------------------------------------------------------------------------

k1 = feval(RHS, u,            xin         );
k2 = feval(RHS, u+dudt*tau/2, xin+k1*tau/2);
xout = xin + tau * k2;
end
