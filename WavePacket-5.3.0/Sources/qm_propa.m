%------------------------------------------------------------------------------
%
% Solves the time-dependent Schroedinger equation to propagate a wave function
% in time.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_propa()

global time;

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('----------------------------------------------------------------------')
util.disp (' Numerical solution of coupled time-dependent Schroedinger            ')
util.disp (' equations using atomic units, i.e. hbar = m_e = e = 1                ')
util.disp (' for quantum systems interacting with electrical fields               ')
util.disp (' (using semiclassical dipole approximation)                           ')
util.disp (' using partial differential equations (pde) solvers                   ')
util.disp (' https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_propa  ')
util.disp ('                                                                      ')
util.disp ('    d                 (        d    )                                 ')
util.disp (' i -- psi( R, t ) = H ( R, -i --, t ) psi( R, t )                     ')
util.disp ('   dt                 (       dR    )                                 ')
util.disp ('                                                                      ')
util.disp (' with psi(R,t)             Wavefunction      (vector, complex)        ')
util.disp (' with H = T + V -iW -F*\mu -F^2/2*\alpha     Hamilton operator        ')
util.disp (' with T = T(-i d/dR)       Kinetic energy    (scalar)                 ')
util.disp (' with V = V(R)             Potential energy  (matrix)                 ')
util.disp (' with W = W(R)             Negative imaginary potential (scalar)      ')
util.disp (' with F = F(t)             Electrical field  (scalar, time-dependent) ')
util.disp (' with \mu = \mu(R)         Dipole moment     (matrix)                 ')
util.disp (' with \alpha = \alpha(R)   Polarizability    (vector)                 ')
util.disp ('                                                                      ')
util.disp ('----------------------------------------------------------------------')
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

%% Beginning of main loop over time steps (step=1 is the initial step )
for step = 1 : time.main.n
    
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
