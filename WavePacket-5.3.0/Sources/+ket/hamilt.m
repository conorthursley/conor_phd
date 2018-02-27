%---------------------------------------------------------------------
%
% Hamiltonian operator acting on wavefunction
%
% Arguments e_x, e_y are components of external electric field
% or half its envelope (if dressed states are used).
% Argument  norm==1 to use normalized Hamiltonian (see Chebychev)
% 
% Original wavefunction is assumed in field psi.dvr.grid_ND (unchanged!)
% Resulting wavefunction is stored in field psi.dvr.new_ND
%
%---------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011      Ulf Lorenz
%
% see the README file for license details.

function hamilt (e_x, e_y, norm)
global psi hamilt space


%% First of all, apply external kinetic energy operators
kinetic_grid = cell(hamilt.coupling.n_eqs, 1);
for m = 1:hamilt.coupling.n_eqs
    kinetic_grid{m} = zeros(size(psi.dvr.grid_ND{1}));
end

if isfield(hamilt, 'kin')
    for n = 1:length(hamilt.kin)
        hamilt.kin{n} = apply(hamilt.kin{n}, false);
        for m = 1:hamilt.coupling.n_eqs
            kinetic_grid{m} = kinetic_grid{m} + psi.dvr.new_ND{m};
        end
    end
end
    
% Apply the kinetic energy operator for each degree of freedom.
% Since the kinetic energy operators are summed up, we have to do
% so as well.
for k = 1:space.size.n_dim
    space.dof{k} = kinetic(space.dof{k}, false);
    for m = 1:hamilt.coupling.n_eqs
        kinetic_grid{m} = kinetic_grid{m} + psi.dvr.new_ND{m};
    end
end
psi.dvr.new_ND = kinetic_grid;
    

% Loop over uncoupled wavefunctions
for m = 1:hamilt.coupling.n_eqs

    % Potential energy (diagonal)
    % if ~isempty ( hamilt.pot.grid_ND{m,m} )
    psi.dvr.new_ND{m} = psi.dvr.new_ND{m} + ...
        hamilt.pot.grid_ND{m,m} .* psi.dvr.grid_ND{m};

    % Permanent dipole moment along x
    if abs(e_x)>0 && ~isempty ( hamilt.d_x.grid_ND{m,m} )
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - ...
            e_x * hamilt.d_x.grid_ND{m,m} .* psi.dvr.grid_ND{m};
    end
    
    % Induced dipole moment along x
    if abs(e_x)>0 && ~isempty ( hamilt.p_x.grid_ND{m,m} )
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - ...
            e_x^2/2 * hamilt.p_x.grid_ND{m,m} .* psi.dvr.grid_ND{m};
    end
    
    % Permanent dipole moment along y
    if abs(e_y)>0 && ~isempty ( hamilt.d_y.grid_ND{m,m} )
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - ...
            e_y * hamilt.d_y.grid_ND{m,m} .* psi.dvr.grid_ND{m}; 
    end
    
    % Induced dipole moment along y
    if abs(e_y)>0 && ~isempty ( hamilt.p_y.grid_ND{m,m} )
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - ...
            e_y^2/2 * hamilt.p_y.grid_ND{m,m} .* psi.dvr.grid_ND{m}; 
    end
    
end
    
% Offdiagonal (coupling) elements of Hamiltonian ( m < n )
% Potential, dipoles but no polarizabilities yet
for m = 1:hamilt.coupling.n_eqs-1
    for n = m+1:hamilt.coupling.n_eqs
        
        % Diabatic potential coupling
       if ~isempty ( hamilt.pot.grid_ND{m,n} )
            psi.dvr.new_ND{m} = psi.dvr.new_ND{m} +      hamilt.pot.grid_ND{m,n}  .* psi.dvr.grid_ND{n};
            psi.dvr.new_ND{n} = psi.dvr.new_ND{n} + conj(hamilt.pot.grid_ND{m,n}) .* psi.dvr.grid_ND{m};
        end
        
        % Transition dipole moment along x
        if abs(e_x)>0 && ~isempty ( hamilt.d_x.grid_ND{m,n} )
            psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - e_x *      hamilt.d_x.grid_ND{m,n}  .* psi.dvr.grid_ND{n};
            psi.dvr.new_ND{n} = psi.dvr.new_ND{n} - conj(e_x * hamilt.d_x.grid_ND{m,n}) .* psi.dvr.grid_ND{m};
        end
        
        % Transition dipole moment along y
        if abs(e_y)>0 && ~isempty ( hamilt.d_y.grid_ND{m,n} )
            psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - e_y *      hamilt.d_y.grid_ND{m,n}  .* psi.dvr.grid_ND{n};
            psi.dvr.new_ND{n} = psi.dvr.new_ND{n} - conj(e_y * hamilt.d_y.grid_ND{m,n}) .* psi.dvr.grid_ND{m};
        end
        
    end
end

% If desired, use normalized Hamiltonian (for Chebychev propagators only)
if (norm==1)
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} - psi.dvr.grid_ND{m} * ( hamilt.range.delta/2 + hamilt.range.min );
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} * 2 / hamilt.range.delta;
    end
end
