%------------------------------------------------------------------------------
%
% Initialize expectation values and their uncertainties
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function expect
global expect space time uncert hamilt

%% Basic setup

% Number of spatial dimensions, additional multiplicative operators
k = space.size.n_dim;
if isfield(space, 'amo')
    p = length(space.amo);
end
    
% Minimum population to output expectation values
if ~isfield(expect, 'min_pop')
    expect.min_pop = 10^-3;
end

%% Expectation values and uncertainties of individual wavefunctions
for m=1:hamilt.coupling.n_eqs

    expect.ind.pop {m} = zeros ( time.main.n, 1 ); % Population

    if isfield(space, 'amo')
        expect.ind.amo {m} = zeros ( time.main.n, p ); % Projection function
    end
    
    expect.ind.dvr {m} = zeros ( time.main.n, k ); % Position
    expect.ind.dvr2{m} = zeros ( time.main.n, k ); % Position^2
    uncert.ind.dvr {m} = zeros ( time.main.n, k ); % Uncertainty
    
    expect.ind.fbr {m} = zeros ( time.main.n, k ); % Momentum
    expect.ind.fbr2{m} = zeros ( time.main.n, k ); % Momentum^2
    uncert.ind.fbr {m} = zeros ( time.main.n, k ); % Uncertainty

    % only used for reduced density plots
    expect.ind.coh{m}  = zeros ( time.main.n, k ); % coherences

    expect.ind.pot {m} = zeros ( time.main.n, 1 ); % Potential energy
    expect.ind.pot2{m} = zeros ( time.main.n, 1 ); % Potential energy^2
    uncert.ind.pot {m} = zeros ( time.main.n, 1 ); % Uncertainty

    expect.ind.kin {m} = zeros ( time.main.n, 1 ); % Kinetic energy
    expect.ind.kin2{m} = zeros ( time.main.n, 1 ); % Kinetic energy^2
    uncert.ind.kin {m} = zeros ( time.main.n, 1 ); % Uncertainty

    expect.ind.all {m} = zeros ( time.main.n, 1 ); % Sum of all energies

end
    
%% Expectation values of total wavefunction

expect.tot.pop = zeros ( time.main.n, 1 ); % Population

if isfield(space, 'amo')
    expect.tot.amo = zeros ( time.main.n, p ); % Projection function
end
    
expect.tot.dvr = zeros ( time.main.n, k ); % Position
expect.tot.fbr = zeros ( time.main.n, k ); % Momentum

expect.tot.coh = zeros ( time.main.n, k ); % coherence

expect.tot.pot = zeros ( time.main.n, 1 ); % Potential energy
expect.tot.kin = zeros ( time.main.n, 1 ); % Kinetic energy
expect.tot.all = zeros ( time.main.n, 1 ); % Sum of pot. and kin. energy

expect.tot.tot = zeros ( time.main.n, 1 ); % Total energy
