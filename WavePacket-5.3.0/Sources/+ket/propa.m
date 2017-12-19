%--------------------------------------------------------------------------
%
% Propagate wavefunction subject to a given Hamiltonian
% using a user-specified PDE solvers for TDSEs
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008,2010 Ulf Lorenz
%
% see the README file for license details.

function propa (step)
global psi hamilt time

% Temporal grids based on sub steps: time, autocorrelation, electric fields
if step==1
    time.sub.short      = time.main.start*time.main.delta;
    time.acf.short      = 1;
    if time.efield.n_pulse>0
        time.efield.short.x = 0;
        time.efield.short.y = 0;
    end
else
    time.sub.short      = [ 0:time.sub.n ]'* time.sub.delta + ...
        (time.main.start+step-2)*time.main.delta;
    time.acf.short      = zeros(time.sub.n, 1);
    if time.efield.n_pulse>0 
        time.efield.short.x = zeros(time.sub.n+1, 1);
        time.efield.short.y = zeros(time.sub.n+1, 1);
    end
end

% Evaluate electric field
if time.efield.n_pulse>0
   time.efield.short = util.efield (time.sub.short); 
end

% Use a propagator of your choice
feval (  time.propa.handle, step );
 
% Apply negative imaginary potential (absorbing boundary conditions)
% using a Trotter splitting
if ~isempty (hamilt.nip.grid_ND)
    for m=1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} .* exp( - hamilt.nip.grid_ND ) ;
    end
end

% Initialize/append temporal grid, electric field along x/y, autocorrelation
if step==1
	time.acf.grid = zeros((time.main.n-1)*time.sub.n+1, 1);
    time.acf.grid(1) = time.acf.short;
    if time.efield.n_pulse>0
		time.efield.grid.x = zeros((time.main.n-1)*time.sub.n+1, 1);
		time.efield.grid.y = zeros((time.main.n-1)*time.sub.n+1, 1);
        time.efield.grid.x(1) = time.efield.short.x;
        time.efield.grid.y(1) = time.efield.short.y;
    end
else
    indexes = 1 + (step-2) * time.sub.n + [1:time.sub.n];
    time.acf.grid(indexes) = time.acf.short;
    if time.efield.n_pulse>0
        time.efield.grid.x(indexes) = time.efield.short.x(2:end);
        time.efield.grid.y(indexes) = time.efield.short.y(2:end);
    end
end

