%------------------------------------------------------------------------------
%
% Solves the time-independent Schroedinger equation to obtain the bound states
% of a system.
%
%------------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_bound()

% Main variables are global throughout;
global time

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('----------------------------------------------------------------------')
util.disp (' Numerical solution of (one or several coupled!) time-independent     ')
util.disp (' Schrödinger equation(s) for bound state problems using Fourier-Grid- ')
util.disp (' Hamiltonian  or other discrete variable representation (DVR) schemes ')
util.disp (' Atomic units (with hbar = m_e = e = 1) are used throughout this code ')
util.disp (' https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_bound  ')
util.disp ('                                                                      ')
util.disp ('   (       d  )                                                       ')
util.disp (' H ( R, -i -- ) psi_n ( R ) = E_n psi_n ( R )                         ')
util.disp ('   (       dR )                                                       ')
util.disp ('                                                                      ')
util.disp (' with psi_n(R)         Wavefunction (n-th eigenvector of H, vector)   ')
util.disp (' with E_n              Eigenenergy (eigenvalue of H)                  ')
util.disp (' with H = T + V -iW    Hamilton operator                              ')
util.disp (' with T = T(R,-i d/dR) Kinetic energy (scalar)                        ')
util.disp (' with V = V(R)         Potential energy (matrix)                      ')
util.disp (' with W = W(R)         Negative imaginary potential (scalar)          ')
util.disp ('                                                                      ')
util.disp ('----------------------------------------------------------------------')
util.disp (' ')

%% Initialize coupling scheme and spatial discretization
init.grid;

% Initialize Hamiltonian operator 
init.hamilt;

% Solve eigenproblem using a finite basis representation of Hamiltonian
init.eigen;

% Initialize expectation values and uncertainties
init.expect;

%% Beginning of main loop over bound states
for step = 1:time.main.n
    
    % Extract psi from columns(!) of eigenvector matrix and normalize
    ket.eigen ( step );
    
    % Expectation values and uncertainties
    ket.expect ( step );
    
    % Logging and plot title
    util.logging ( step );       
             
    % Plot densities and expectation values
    plot.wigner ( step );
    plot.psi ( step );

    % save the wave function
    ket.save ( step );
    
% End of main loop     
end

% Output clock/date/time
util.clock;

end
