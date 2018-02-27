%------------------------------------------------------------------------------
%
% Restarts an interrupted calculation at a given time step.
%
% If you run a calculation and save the results, and if this calculation is
% interrupted unexpectedly, you can use this script to restart the calculation.
% After you have regenerated some global information using the regenerate()
% script, this script starts the calculation from some time step that was
% already saved.
%
% Input parameters are the directory where the saved calculations reside
% (ket.save.dir of the interrupted calculation), the file name of the saved files
% (ket.save.file), and the step from which you restart the calculation.
%
% Note that a handful of variables will be broken, since we leave out all steps
% in-between. These are all of expect and uncertain (they have gaps where we
% skipped the calculation of expectation values), time.acf.grid and
% time.efield.grid.* (wrong indices), time.freq.grid and time.spec.grid
% (inconsistent data)
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2012 Ulf Lorenz
%
% see the README file for license details.

function qm_rerun(dir, file, start_step)

global time;

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' Numerical solution of coupled time-dependent Schr�dinger    ')
util.disp (' equations using atomic units, i.e. hbar = m_e = e = 1       ')
util.disp (' for quantum systems interacting with electrical fields      ')
util.disp (' (using semiclassical dipole approximation)                  ')
util.disp (' using partial differential equations (pde) solvers          ')
util.disp ('                                                             ')
util.disp ('    d                 (        d    )                        ')
util.disp (' i -- psi( R, t ) = H ( R, -i --, t ) psi( R, t )            ')
util.disp ('   dt                 (       dR    )                        ')
util.disp ('                                                             ')
util.disp (' with psi(R,t)           Wavefunction      (vector, complex) ')
util.disp (' with H = T + V -iW -F*D Hamilton operator                   ')
util.disp (' with T = T(-i d/dR)     Kinetic energy    (scalar)          ')
util.disp (' with V = V(R)           Potential energy  (matrix)          ')
util.disp (' with W = W(R)           Negative imag. potential (scalar)   ')
util.disp (' with D = D(R)           Dipole moment     (matrix)          ')
util.disp (' with F = F(t)           Electrical field  (scalar, time-dep)')
util.disp ('                                                             ')
util.disp ('-------------------------------------------------------------')
util.disp (' ')


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

%% Run one step for the setup
for step = 1
    
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
    
% End of main loop
end


% Load the wave function
context = ket.load(dir, file);
ket.load(context, start_step, true);
clear context


% And continue from where we left off.
for step = start_step+1 : time.main.n
    
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
    
% End of main loop
end

% Output clock/date/time
util.clock;

end