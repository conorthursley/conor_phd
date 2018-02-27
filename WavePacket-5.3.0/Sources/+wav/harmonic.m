%------------------------------------------------------------------------------
%
% This function creates the initial state 
% as an eigenstate of a Harmonic oscillator.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% see the README file for license details.

function init_grid = harmonic(dir)
global space psi

%% Set default values and check
if ~isfield(psi.dof{dir}, 'm_r')
        psi.dof{dir}.m_r = space.dof{dir}.mass;
end
if ~isfield(psi.dof{dir}, 'r_e')
        psi.dof{dir}.r_e = 0;
end
if ~isfield(psi.dof{dir},'n_q')
    psi.dof{dir}.n_q = 0;
end
if (psi.dof{dir}.n_q<0)
    util.error ('quantum number must not be negative')
end


% Use either v_2 or omega from the input file; calculate missing quantity
if isfield(psi.dof{dir}, 'v_2') && ~isfield(psi.dof{dir}, 'omega')
    psi.dof{dir}.omega = sqrt(psi.dof{dir}.v_2 / psi.dof{dir}.m_r);
end
if isfield(psi.dof{dir}, 'omega') && ~isfield(psi.dof{dir}, 'v_2')
    psi.dof{dir}.v_2 = psi.dof{dir}.omega^2 * psi.dof{dir}.m_r;
end

%% Some output
util.disp (' ')
util.disp ('*******************************************************')
util.disp ( ['Initial wavefunction for DOF :' int2str(dir)] )
util.disp ('   ' )
util.disp ('Harmonic oscillator eigenstate')
util.disp ( ['Angular Frequency            : ' num2str(psi.dof{dir}.omega)] )
util.disp ( ['corresp. Force constant      : ' num2str(psi.dof{dir}.v_2)] )
util.disp ( ['(Reduced) mass               : ' num2str(psi.dof{dir}.m_r)] )
util.disp ( ['Equilibrium position         : ' num2str(psi.dof{dir}.r_e)] )
util.disp ( ['Quantum number (eigenstate)  : ' int2str(psi.dof{dir}.n_q)] )


%% Set up the grid

% HO eigenstate \Psi(x) = exp(-m\omega/2 * x^2) * H_n(\sqrt(m\omega) * x)
factor    = psi.dof{dir}.m_r * psi.dof{dir}.omega;
position  = space.dvr.grid_ND{dir} - psi.dof{dir}.r_e;
init_grid = exp(- factor/2 * position.^2) ...
            .* math.hermite(sqrt(factor) * position, psi.dof{dir}.n_q);
