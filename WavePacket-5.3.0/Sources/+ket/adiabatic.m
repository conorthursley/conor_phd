%-------------------------------------------------------------------------
%
% Transform (matrix-valued) potential and (vector-valued) wavefunction
% from diabatic to adiabatic representation or back. It is important 
% that the dia2adi transformation has to be called BEFORE adi2dia
% because the latter uses data stored by the former.
%
% Note that within WavePacket the adiabatc picture is used only for the 
% purpose of visualizing wavefunctions and calculating corresponding 
% expectation values but NOT for the propagation itself! This is
% because of the problems associated with the (near) singularities 
% of the kinetic (derivative) couplings at (avoided) crossings and 
% intersections.
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
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function adiabatic ( direction )
global hamilt psi space time

% Use this routine only if adiabatic representation is desired
if strcmpi(hamilt.coupling.representation,'dia')
    return
end

% Instantaneous electric field present?
if ~isfield(time, 'efield') || ~isfield(time.efield, 'n_pulse') || time.efield.n_pulse==0
    e_x = 0;
    e_y = 0;
else
    e_x = time.efield.short.x(1);
    e_y = time.efield.short.y(1);
end


%% Diabatic to adiabatic ("forward transformation")
if strcmpi(direction,'dia2adi')

    % Save diabatic wave functions and potential  matrices (cell arrays!)
    psi.dvr.dia    = psi.dvr.grid_ND;
    hamilt.pot.dia = hamilt.pot.grid_ND;

    switch hamilt.coupling.n_eqs

       % Light-dressed potential
       case 1

            dia = hamilt.pot.dia{1,1};
            
            % Permanent dipole moment along x
            if abs(e_x)>0 && ~isempty ( hamilt.d_x.grid_ND{1,1} )
                dia = dia - e_x * hamilt.d_x.grid_ND{1,1};
            end
            
            % Induced dipole moment along x
            if abs(e_x)>0 && ~isempty ( hamilt.p_x.grid_ND{1,1} )
                dia = dia - e_x^2/2 * hamilt.p_x.grid_ND{1,1};
            end
            
            % Permanent dipole moment along y
            if abs(e_y)>0 && ~isempty ( hamilt.d_y.grid_ND{1,1} )
                dia = dia - e_y * hamilt.d_y.grid_ND{1,1};
            end
            
            % Induced dipole moment along y
            if abs(e_y)>0 && ~isempty ( hamilt.p_y.grid_ND{1,1} )
                dia = dia - e_y^2/2 * hamilt.p_y.grid_ND{1,1};
            end
            
            % "Adiabatic" potential
            hamilt.pot.grid_ND{1,1} = dia;
            
            % Wavefunctions unchanged
            psi.dvr.grid_ND{1} = psi.dvr.dia{1};
        
        % Analytic diagonalization for two coupled Schrödinger equations
        case 2

            % Diabatic potentials
            dia1 = hamilt.pot.dia{1,1};
            dia2 = hamilt.pot.dia{2,2};

            % Permanent dipole moments along x
            if abs(e_x)>0
                if ~isempty ( hamilt.d_x.grid_ND{1,1} )
                    dia1 = dia1 - e_x * hamilt.d_x.grid_ND{1,1};
                end
                if ~isempty ( hamilt.d_x.grid_ND{2,2} )
                    dia2 = dia2 - e_x * hamilt.d_x.grid_ND{2,2};
                end
            end

            % Induced dipole moments along x
            if abs(e_x)>0
                if ~isempty ( hamilt.p_x.grid_ND{1,1} )
                    dia1 = dia1 - e_x^2/2 * hamilt.p_x.grid_ND{1,1};
                end
                if ~isempty ( hamilt.p_x.grid_ND{2,2} )
                    dia2 = dia2 - e_x^2/2 * hamilt.p_x.grid_ND{2,2};
                end
            end

            % Permanent dipole moments along y
            if abs(e_y)>0
                if ~isempty ( hamilt.d_y.grid_ND{1,1} )
                    dia1 = dia1 - e_y * hamilt.d_y.grid_ND{1,1};
                end
                if ~isempty ( hamilt.d_y.grid_ND{2,2} )
                    dia2 = dia2 - e_y * hamilt.d_y.grid_ND{2,2};
                end
            end

            % Induced dipole moments along y
            if abs(e_y)>0
                if ~isempty ( hamilt.p_y.grid_ND{1,1} )
                    dia1 = dia1 - e_y^2/2 * hamilt.p_y.grid_ND{1,1};
                end
                if ~isempty ( hamilt.p_y.grid_ND{2,2} )
                    dia2 = dia2 - e_y^2/2 * hamilt.p_y.grid_ND{2,2};
                end
            end

            % Sum, difference of diabatic potentials
            eta   = (dia1 + dia2)/2; % half trace
            delta = (dia1 - dia2)/2; % half gap

            % Three possible coupling mechanisms
            beta  =  zeros(size(space.dvr.grid_ND{1}));

            % Diabatic potential coupling
            if ~isempty(hamilt.pot.dia{1,2})
                beta  =  beta + hamilt.pot.dia{1,2};
            end

            % Transition dipole moments along x
            if abs(e_x)>0 && ~isempty(hamilt.d_x.grid_ND{1,2})
                beta = beta - e_x * hamilt.d_x.grid_ND{1,2};
            end

            % Transition dipole moments along y
            if abs(e_y)>0 && ~isempty(hamilt.d_y.grid_ND{1,2})
                beta = beta - e_y * hamilt.d_y.grid_ND{1,2};
            end

            % Polar coordinates: "Radius", mixing angle
            theta = atan2 ( beta , delta );
            rho   = sqrt ( delta.^2 + beta.^2 );

            % Adiabatic potential energy surfaces: 1=lower, 2=upper
            hamilt.pot.grid_ND{1,1} = eta - rho;
            hamilt.pot.grid_ND{2,2} = eta + rho;
            hamilt.pot.grid_ND{1,2} = zeros(size(hamilt.pot.grid_ND{1,2}));

            % Adiabatic wavefunctions
            psi.dvr.grid_ND{1} = exp(1i*theta/2) .* ( - sin(theta/2).*psi.dvr.dia{1} + cos(theta/2).*psi.dvr.dia{2} );
            psi.dvr.grid_ND{2} = exp(1i*theta/2) .* ( + cos(theta/2).*psi.dvr.dia{1} + sin(theta/2).*psi.dvr.dia{2} );

        % Numerical diagonalization for more than two coupled Schrödinger equations
        otherwise

            psi.vec = zeros(hamilt.coupling.n_eqs,1);
            for gr=1:space.size.n_tot
                pot_mat = zeros(hamilt.coupling.n_eqs);

                % Extract diabatic wavefunctions, diabatic potentials from cell arrays
                for m=1:hamilt.coupling.n_eqs
                    psi.vec(m)   = psi.dvr.dia{m}(gr);
                    pot_mat(m,m) = hamilt.pot.dia{m,m}(gr);

                    % Permanent dipole moments along x
                    if abs(e_x)>0 && ~isempty ( hamilt.d_x.grid_ND{m,m} )
                        pot_mat(m,m) = pot_mat(m,m) - e_x * hamilt.d_x.grid_ND{m,m}(gr);
                    end

                    % Induced dipole moments along x
                    if abs(e_x)>0 && ~isempty ( hamilt.p_x.grid_ND{m,m} )
                        pot_mat(m,m) = pot_mat(m,m) - e_x^2/2 * hamilt.p_x.grid_ND{m,m}(gr);
                    end

                    % Permanent dipole moments along y
                    if abs(e_y)>0 && ~isempty ( hamilt.d_y.grid_ND{m,m} )
                        pot_mat(m,m) = pot_mat(m,m) - e_y * hamilt.d_y.grid_ND{m,m}(gr);
                    end

                    % Induced dipole moments along y
                    if abs(e_y)>0 && ~isempty ( hamilt.p_y.grid_ND{m,m} )
                        pot_mat(m,m) = pot_mat(m,m) - e_y^2/2 * hamilt.p_y.grid_ND{m,m}(gr);
                    end


                    % Three possible coupling mechanisms
                    for n=m+1:hamilt.coupling.n_eqs
                        % Diabatic potential coupling
                        if ~isempty(hamilt.pot.dia{m,n})
                            pot_mat(m,n) = pot_mat(m,n) + hamilt.pot.dia{m,n}(gr);
                        end

                        % Transition dipole moments along x
                        if abs(e_x)>0 && ~isempty(hamilt.d_x.grid_ND{m,n})
                            pot_mat(m,n) = pot_mat(m,n) - e_x * hamilt.d_x.grid_ND{m,n}(gr);
                        end

                        % Transition dipole moments along y
                        if abs(e_y)>0 && ~isempty(hamilt.d_y.grid_ND{m,n})
                            pot_mat(m,n) = pot_mat(m,n) - e_y * hamilt.d_y.grid_ND{m,n}(gr);
                        end
                        
                        % Make potential matrix Hermitian 
                        pot_mat(n,m) = conj(pot_mat(m,n)); 
                    end
                end

                % Transform wavefunctions, potentials to adiabatic picture
                [eig_vec, eig_val] = eig(pot_mat);
                psi.vec = eig_vec' * psi.vec;

                % Save adiabatic wavefunctions, potentials in cell arrays
                for n=1:hamilt.coupling.n_eqs
                    psi.dvr.grid_ND   {n}  (gr) = psi.vec(n);
                    hamilt.pot.grid_ND{n,n}(gr) = eig_val(n,n);
                end


            end
    end

    %% Adiabatic to diabatic ("back transformation")
elseif strcmpi(direction,'adi2dia')

    % Retrieve diabatic wavefunctions and potential matrices (cell arrays!)
    psi.dvr.grid_ND    = psi.dvr.dia;
    hamilt.pot.grid_ND = hamilt.pot.dia;

end
