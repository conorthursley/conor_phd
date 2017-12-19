function qm_init
global hamilt psi space time control 

util.disp ( '***********************************' )
util.disp ( 'Optimally driven double-well system' )
util.disp ( '                                   ' )
util.disp ( 'see J. Werschnik and E. K. U. Gross' )
util.disp ( 'J. Phys. B 40, R175  (2007)        ' )
util.disp ( 'doi:10.1088/0953-4075/40/18/R01    ' )
util.disp ( '***********************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max = +10.0;              % Upper bound of grid

% Hamiltonian operator 
hamilt.pot.handle = @pot.taylor;         % Quartic double well potential
hamilt.pot.v{1,1}  = [0; -2/4; +6/256; +24/64]; 

hamilt.dip.handle = @dip.taylor;         % Dipole moment: Taylor series
hamilt.dip.x.v{1,1} = 1;                 % Linear dipole moment, slope 1

hamilt.eigen.cutoff = 0.0;               % Cut-off entries of Hamiltonian matrix
hamilt.truncate.min = -05.0;             % Lower truncation of energy
hamilt.truncate.max = +20.0;             % Upper truncation of energy

% Calculate and save (bound) eigen states (==> qm_bound)
psi.eigen.start = 0;                     % Lower index
psi.eigen.stop  = 20;                    % Upper index
psi.save.export = true;                  % Toggle export
psi.save.dir = pwd;                      % Directory/folder
psi.save.file = 'bound';                 % File name (.mat)

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 040;                   % Index of final time step
time.main.delta = 010;                   % Size of main time steps 
time.sub.n      = 050;                   % Number of substeps per main step

% Electric field pulses: Initial guess ==> qm_optimal
% Works best with ampli<0 (otherwise fluence too high)
time.efield.dressed = false;
time.efield.shape   = 'recta';           % Shape of envelope
time.efield.polar   = 0;                 % Polarization angle [rad]
time.efield.delay   = 200;               % Time delay of pulse center
time.efield.fwhm    = 400;               % Pulse length
time.efield.ampli   = -0.05;             % Amplitude of electric field         
time.efield.frequ   = 0;                 % Carrier frequency
time.efield.phase   = 0;                 % Phase

% frequency resolved optical gating
time.frog.choice    = 'Gauss';          % choice of methods
time.frog.zoom      = 15;                % zoom factor for frequency axis
time.frog.width     = 20;                % zooming frequency axis

% Initial density
control.initial.choice = 'pure';         % starting from a pure state
control.initial.pure = 0;                % selecting ground state

% Select ODE solvers, parameters
control.solvers.handle1 = @ode.RuKu6;    % Runge-Kutta or other integrators
control.solvers.handle2 = @ode45;        % Runge-Kutta: Dormand-Prince (built-in)
control.solvers.reltol = 1e-6;           % Relative error tolerance of solution

% Define output observables and choose control targets
control.observe.types = 'ovl';           % choose types of observables
control.observe.choices={0 1 2:20}; 
control.observe.choices={0 1 2:20}; 
control.observe.labels={'Left well' 'Right well' 'Delocalized'};
control.observe.targets=1:3; 

% Optimal control theory
control.optimal.terminal = 2;            % Which observable to be optimized
control.optimal.tolerance = 1e-4;        % Threshold terminating iteration
control.optimal.max_iter = 050;          % Max. number of iteration steps
control.optimal.alpha = 2.20;            % Penalty factor for laser fluence
control.optimal.eta  = 0.90;             % Calculate optimal backward field
control.optimal.zeta = 0.10;             % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields

% Plot settings
control.plot.uxy = true;         % Plot u(t), x(t), y(t)
control.plot.j12 = true;         % Plot j_1(t), j_2(t), and total
control.plot.psd = true;         % Plot power spectral density
control.plot.mov = true;         % Animation of u(t), x(t), y(t)

