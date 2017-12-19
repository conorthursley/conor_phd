%-------------------------------------------------------------------------
%
% Transform (vector-valued) wavefunction from adiabatic to 
% diabatic representation. Within WavePacket, this routine 
% is used _only_ for the purpose of transforming an initial 
% wavefunction in an adiabatic representation to the diabatic
% representation for the first time. If you wish to transform
% the wave function later, use ket.adiabatic.
%
% Note that this function has not yet been implemented for more
% than two Schroedinger equations.
%
% For the special case of 2 coupled Schrödinger equations, 
% the transformations can be calculated analytically 
% 
% Diabatic: V = ( alpha  beta )
%               ( beta  gamma )
%
% delta = (alfa-gamma)/2
%   eta = (alfa+gamma)/2
%   rho = sqrt(delta^2+beta^2)>0
% theta = atan ( beta/delta )
%
% Adiabatic: E = (eta-rho    0    )
%                (   0    eta+rho )
%
% dia2adi: S = ( -sin(theta/2)  +cos(theta/2) )
%              ( +cos(theta/2)  +sin(theta/2) )
%
% The first and second column of S contain the eigenvectors
% corresponding to the lower (label 1) and upper (label 2)
% eigenvalue contained in the diagonal entries of matrix E.
% Note that in this case we have S^-1 = S^+ = S, i.e. the
% 'dia2adi' and 'adi2dia' transformation are described by 
% the same matrix.
%
% Obviously, this transformation results in double-valued
% adiabatic states upon increasing theta from 0 to 2*pi, e.g.
% upon encircling a conical intersection of two adiabatic
% potential energy surfaces. This can be remidied by intro-
% ducing complex-valued adiabatic states, i.e. by multiplying
% S with exp(i*phi) or multiplying S^+ with exp(-i*phi). The
% geometric phase has to be chosen as phi = n*theta/2 where
% n is an odd integer. Note that this is equivalent to the 
% introduction of an appropriately chosen vector potential 
% governing the dynamics of wavepackets along the adiabatic 
% potential energy surfaces. (Alternatively, one could compensate
% the double-valuedness of the adiabatic states by enforcing 
% doubled-valued boundary conditions for the wavefunctions, too.)
%
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

function diabatic
global hamilt psi space

% Check number of (coupled) equations
switch hamilt.coupling.n_eqs
    case 1
        return

    case 2
        % Difference, coupling, mixing angle
        delta = (hamilt.pot.grid_ND{1,1} - hamilt.pot.grid_ND{2,2})/2;
        beta  =  hamilt.pot.grid_ND{1,2};
        theta = atan2 ( beta , delta );

        % Diabatic wavefunctions
        adi = psi.dvr.grid_ND;
        psi.dvr.grid_ND{1} = exp(-i*theta/2) .* ( - sin(theta/2).*adi{1} + cos(theta/2).*adi{2} );
        psi.dvr.grid_ND{2} = exp(-i*theta/2) .* ( + cos(theta/2).*adi{1} + sin(theta/2).*adi{2} );
            
    otherwise
        % Basically, cut&paste from ket.adiabatic.
        % For each grid point the electronic problem, but transform the wavefunction
        % "backwards".
        psi.adi = psi.dvr.grid_ND;
        psi.vec = zeros(hamilt.coupling.n_eqs,1);

        for gr=1:space.size.n_tot
            pot_mat = zeros(hamilt.coupling.n_eqs);

            for m=1:hamilt.coupling.n_eqs
                psi.vec(m)   = psi.adi{m}(gr);
                pot_mat(m,m) = hamilt.pot.grid_ND{m,m}(gr);

                for n=m+1:hamilt.coupling.n_eqs
                    if ~isempty(hamilt.pot.grid_ND{m,n})
                        pot_mat(m,n) = hamilt.pot.grid_ND{m,n}(gr);
                        pot_mat(n,m) = conj(pot_mat(m,n)); 
                    end
                end
            end

            [eig_vec, eig_val] = eig(pot_mat);
            psi.vec = eig_vec * psi.vec;

            % Write the result back
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr.grid_ND{m}(gr) = psi.vec(m);
            end
        end
end
