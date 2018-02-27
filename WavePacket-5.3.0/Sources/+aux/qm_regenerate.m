%------------------------------------------------------------------------------
%
% Regenerates some global data needed to restart a saved calculation.
%
% The problem: If you save a calculation and this calculation is interrupted,
% e.g., by a crash, then all the global variables are not written out. They are
% needed, however, to load the wave functions again (ket.load() is written in
% such a way). This function regenerates the global variables, so that the
% calculation can later be started by the rerun() function.
%
% To use this function, setup the calculation, but then instead of qm_propa() et al.
% call regenerate().
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% see the README file for license details.

function qm_regenerate()

global time;

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

%% Initialize coupling scheme and spatial/temporal discretization
init.grid;

% Initialize the electric field
init.efield;

% Initialize Hamiltonian operator 
init.hamilt;

% Initialize temporal discretization
init.timesteps;

% Initialize wave function
init.psi;

% Initialize expectation values and uncertainties
init.expect;

%% Beginning of main loop over time steps (step=1 is the initial step )
step = 1;
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    ket.propa ( step );

    % Transform to adiabatic representation (if desired)
    ket.adiabatic ( 'dia2adi' );
    
    % Expectation values and uncertainties
    ket.expect ( step );
    
    % Logging and plot title
    util.logging ( step );       

    % Get spectrum as Fourier transform of autocorrelation
    ket.spectrum ( step );

    % Plot densities and expectation values
    plot.wigner ( step );
    plot.psi ( step );

    % Transform back to diabatic representation
    ket.adiabatic ( 'adi2dia' );

    % Store the wave function.
    ket.save ( step );
    

% This should save the missing global variables
ket.save( time.main.n );


% Output clock/date/time
util.clock;
