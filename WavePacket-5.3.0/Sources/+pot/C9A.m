%------------------------------------------------------------------------------
%
% Potential energy function for the C9A molecule. Since the potentials
% that are fitted from spectroscopic data are completely nuts 
% (the coeffcients by Monte et al. give _something_, but not the
% correct potentials), we use a fourth order potential in the torsion
% angle. The diabatic coupling between the excited state S1 and Sx
% potentials is given by a Gaussian.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function C9A

global hamilt space

%% Error checking and some setup
if space.size.n_dim ~= 1
    util.error ('Potential only valid for a single dimension!')
end

if isa(space.dof{1}, 'grid.legendre')
    % Legendre grid gives angles as cos(phi)
    phigrid = acos(space.dvr.grid_ND{1}) - pi/2;
else
    % other grids have the angle as coordinate
    phigrid = space.dvr.grid_ND{1} - pi/2;
end

% 1,2, 3 refers to the ground (S0), bright (S1) or dark (Sx) state
if ~isfield(hamilt.pot, 'state')
    switch hamilt.coupling.n_eqs
        case 1
            hamilt.pot.state = 1;
        case 2
            hamilt.pot.state = [2 3];
        case 3
            hamilt.pot.state = [1 2 3];
        otherwise
            util.error('This potential is only for up to three surfaces')
    end
end

if length(hamilt.pot.state) ~= hamilt.coupling.n_eqs
    util.error('Number of potential states does not match number of equations!')
end

%% Output
util.disp ('***********************************************')
util.disp ('Potential for 9-(N-carbazolyl) anthracene (C9A)')
util.disp ('S0, S1 and "dark" Sx state. Potentials taken   ')
util.disp ('from Manz et al. Z.Phys.D 34:111               ')
util.disp ('Monte et al. JChemPhys 98:2580                 ')
util.disp ('***********************************************')
util.disp ( ['State(s)    : ' num2str(hamilt.pot.state)])
util.disp ('')


%% Parameter setup.

% For each state, we setup the parameters first. The potentials are
% V(x) = V_0 + V_2/2 (x-pi/2)^2 + V_4/24 (x-pi/2)^4
% V{i}(k) gives the parameter of the i-th surface, k=1,2,3 translates
% to V_0,V_2,V_4. All values in cm^-1, and angles in rad. Case A for the
% Sx potential.
V = {[17 -1430.86 180648.43], [], [25806 19094 0]};

% coupling constants
C   = 22 * 4.5563e-6;       % coupling strength in atomic units
a   = 0.056;                % width of Gaussian
phi0= 0.227;                % center of Gaussian


%% Now calculate the potentials
S1 = [];
Sx = [];

for m = 1:hamilt.coupling.n_eqs
    state = hamilt.pot.state(m);
    if state == 2
        % Here, we use the potential of Monte et al. because a quartic fit gives
        % far too bad results (harmonic frequency is off by a factor of sqrt(2)).
        % On the other hand, the parameters in Monte et al. are completely useless
        % for the groundstate (no idea _what_ they did, though). The parameters
        % for the S1 state also have a wrong sign, but seem to be correct otherwise.
        % The potential is a fit V = \sum_n V_n/2 (1 - cos(n*phi)), of course with
        % a wrong offset that we have to add by hand.
        phi = phigrid + pi/2;
        coeffs = [ 6712 -2777 -6963 -2007 1692 373  -4688 -7518 ...
                  -5501 -602   3348  4539 3511 1873  662   140];
        hamilt.pot.grid_ND{m,m} = zeros(size(phi));

        for ii = 1:length(coeffs)
            hamilt.pot.grid_ND{m,m} = hamilt.pot.grid_ND{m,m} ...
                    + coeffs(ii)/2 * (1 - cos(2*ii*phi));
        end

        hamilt.pot.grid_ND{m,m} = -hamilt.pot.grid_ND{m,m} - 193.3 + 25893.3;
    else
        hamilt.pot.grid_ND{m,m} = V{state}(1) + V{state}(2)/2 * phigrid.^2 ...
                + V{state}(3)/24 * phigrid.^4;
    end

    % convert to atomic units
    hamilt.pot.grid_ND{m,m} = hamilt.pot.grid_ND{m,m} * 4.5563e-6;


    if state == 2
        S1 = m;
    elseif state == 3
        Sx = m;
    end
end

if ~isempty(S1) && ~isempty(Sx)
    phiplus = phigrid + phi0;
    phiminus= phigrid - phi0;

    hamilt.pot.grid_ND{min(S1,Sx), max(S1,Sx)} = C * (...
        exp(-0.5 * (phiplus/a).^2) + exp(-0.5 * (phiminus/a).^2));
end
