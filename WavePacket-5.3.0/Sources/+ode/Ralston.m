function xout = Ralston(RHS,xin,u,dudt,tau)
%--------------------------------------------------------------------------
% Ralston explicit integrator
% https://en.wikipedia.org/???
%--------------------------------------------------------------------------

k1 = feval(RHS, u,              xin);
k2 = feval(RHS, u+dudt*tau*2/3, xin+k1*tau*2/3);

xout = xin + tau * (k1 + 3*k2) / 4;

end
