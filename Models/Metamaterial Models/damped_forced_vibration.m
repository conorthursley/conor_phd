function [amp,phase] = damped_forced_vibration(D,M,f,omega)

% Function to calculate steady state amplitude of
% a forced linear system.
% D is 2nx2n the stiffness/damping matrix
% M is the 2nx2n mass matrix
% f is the 2n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector ‘amp’, giving the amplitude of
% each degree of freedom, and a second vector ‘phase’,
% which gives the phase of each degree of freedom

Y0 = (D+M.*1i*omega)\f;  % The i here is sqrt(-1)
% We dont need to calculate Y0bar - we can just change the sign of
% the imaginary part of Y0 using the 'conj' command
for j =1:length(f)/2
    amp(j) = sqrt(Y0(j)*conj(Y0(j)));
    phase(j) = log(conj(Y0(j))/Y0(j))/(2*i);
end

end