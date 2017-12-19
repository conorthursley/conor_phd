%--------------------------------------------------------------------------
%
% Matrix (eigen) representations of important operators using
% output from previous qm_bound calculation
%
% This function takes as input the path to a save file similar to qm_movie.
% It then assumes that each entry in the save file is an eigenstate of some
% Hamiltonian, and calculates the eigenenergies, matrix elements of dipole
% and/or polarizability along x and/or y (if applicable) and optionally 
% also matrix elements of the system/bath coupling and/or additional multi-
% plicative operators. These vectors/matrices are combined in structure tise 
% and are written to file "tise.mat" in the current working directory. 
%
% An optional third parameter can be supplied to define a cut-off; matrix
% elements are set to zero if their absolute values are below the cut-off, 
% which keeps the output files more readable (and could be used for sparse
% representations later on).
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014-17 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%               
%
% see the README file for license details.

function qm_matrix(savedir, savefile, cutoff)

global control hamilt psi space time

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Provide default values for missing input arguments
if nargin<3
    cutoff=0;
end

if nargin<2
    savefile = 'bound';
end

if nargin<1
    savedir = pwd;
end

util.disp('-------------------------------------------------------------');
util.disp(' Matrix (eigen) representations of important operators       ');
util.disp(' by numerical integration using DVR (quadrature) schemes     ');
util.disp(' https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_matrix/                                                            ');
util.disp('                                                             ');
util.disp([' Reading eigenvectors from file(s): ' savefile '.mat, ' savefile '_0.mat, ...']);
util.disp([' in directory: ' savedir]);
util.disp('-------------------------------------------------------------');
util.disp(' ')

%% Load the input
context = ket.load(savedir, savefile, true);

%% Preallocate fields (ham, dip, pol, sbc, amo) of structure 'tise' (if applicable)
tise.ham = zeros(time.main.n, 1);
util.disp('Energies: \delta_ij E_i = <i|H|j> ')
util.disp(' ')

if isfield(hamilt,'d_x') && ~isempty(hamilt.d_x.grid_ND{1,1})
    tise.d_x = zeros(time.main.n);
    util.disp('Dipole moments: \mu_{x,ij} = <i|\mu_x|j> ')
    util.disp(' ')
end

if isfield(hamilt,'d_y') && ~isempty(hamilt.d_y.grid_ND{1,1})
    hamilt.d_y.grid_ND
    tise.d_y = zeros(time.main.n);
    util.disp('Dipole moments: \mu_{y,ij} = <i|\mu_y|j> ')
    util.disp(' ')
end

if isfield(hamilt,'p_x') && ~isempty(hamilt.p_x.grid_ND{1,1})
    tise.p_x = zeros(time.main.n);
    util.disp('Polarizabilities: \alpha_{x,ij} = <i|\alpha_x|j> ')
    util.disp(' ')
end

if isfield(hamilt,'p_y') && ~isempty(hamilt.p_y.grid_ND{1,1})
    tise.p_y = zeros(time.main.n);
    util.disp('Polarizabilities: \alpha_{y,ij} = <i|\alpha_y|j> ')
    util.disp(' ')
end

if isfield(hamilt,'sbc') && ~isempty(hamilt.sbc.grid_ND{1,1})
    tise.sbc = zeros(time.main.n);
    util.disp('System-bath coupling: \chi_{ij} = <i|\chi|j> ')
    util.disp(' ')
end

if isfield(space, 'amo') && strcmpi(control.observe.types,'amo')
    tise.amo = cell (length(space.amo),1);
    for p = 1:length(space.amo)
        tise.amo{p} = zeros(time.main.n);
        util.disp(['Additional multiplicative operators: O_{p,ij} = <i|O_p|j> with O being: ' space.amo{p}.label ])
        util.disp(' ')
    end
end

%% Calculate all matrix elements
for bra_index = 1:time.main.n
	context = ket.load(context, bra_index, true);
	bra_state = psi.dvr.grid_ND;
    
    % calculate the eigenenergies
    ket.hamilt(0,0,0); % provides H|psi> in psi.dvr.new_ND
    for m = 1:hamilt.coupling.n_eqs
        tise.ham(bra_index) = ket.dot(bra_state, psi.dvr.new_ND);
    end
    
  
	% double loop for matrix elements to be calculated as "sandwiches" 
	for ket_index = 1:time.main.n
		context = ket.load(context, ket_index, true);
		ket_state = psi.dvr.grid_ND;

        % dipole moment(x)
		if isfield (tise, 'd_x')
            tise.d_x(bra_index, ket_index) = ket.sandwich(bra_state, hamilt.d_x.grid_ND, ket_state);
        end

        % dipole moment(y)
		if isfield (tise, 'd_y')
            tise.d_y(bra_index, ket_index) = ket.sandwich(bra_state, hamilt.d_y.grid_ND, ket_state);
        end
        
        % polarization(x)
        if isfield (tise, 'p_x')
            tise.p_x(bra_index, ket_index) = ket.sandwich(bra_state, hamilt.p_x.grid_ND, ket_state);
        end
        
        % polarization(y)
        if isfield (tise, 'p_y')
            tise.p_y(bra_index, ket_index) = ket.sandwich(bra_state, hamilt.p_y.grid_ND, ket_state);
        end
        
        % system-bath coupling
        if isfield (tise, 'sbc')
            tise.sbc(bra_index, ket_index) = ket.sandwich(bra_state, hamilt.sbc.grid_ND, ket_state);
        end
        
        % additional multiplicative operators
        if isfield(tise, 'amo')
            for p = 1:length(space.amo)
                tise.amo{p}(bra_index, ket_index) = ket.sandwich(bra_state, space.amo{p}.grid_ND, ket_state);
            end
        end
        
    end
end

%% Optionally truncate matrix elements
if cutoff>0
	if isfield (tise, 'd_x')
        tise.d_x(abs(tise.d_x) < cutoff) = 0;
    end
	if isfield (tise, 'd_y')
        tise.d_y(abs(tise.d_y) < cutoff) = 0;
    end
	if isfield (tise, 'p_x')
        tise.p_x(abs(tise.p_x) < cutoff) = 0; 
    end
	if isfield (tise, 'p_y')
        tise.p_y(abs(tise.p_y) < cutoff) = 0; 
    end
    if isfield (tise, 'sbc')
        tise.sbc(abs(tise.sbc) < cutoff) = 0; 
    end
    if isfield(tise, 'amo')
        for p = 1:length(space.amo)
            tise.amo{p}(abs(tise.amo{p}) < cutoff) = 0;
        end
    end

end

%% Matrix/vector representations and labels of observables (to be used in qm_abncd)
tise.lab = cell(length(control.observe.choices),1);
tise.obs = control.observe.types;

switch lower(control.observe.types)
    
    % Additional multiplicative operators
    case 'amo'
        if ~isfield (tise, 'amo')
            util.error ('No additional multiplicative operators defined')
        end
        
        % Only observables corresponding to *one* operator
        for len=1:length (control.observe.choices)
            if length(control.observe.choices{len})>1
                util.error('Combinations of more than one observables not *yet* implemented')
            end
        end

        % If not specified otherwise, observables will be labeled like operators
        if ~isfield (control.observe,'labels')
            for len=1:length (control.observe.choices)
                control.observe.labels{len} = space.amo{len}.label;
            end
        end
            
        % Set labels and observables ==> structure "tise"
        tise.mat = cell(length(control.observe.choices),1);
        for len=1:length (control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            util.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' tise.lab{len}])
            tise.mat{len} = tise.amo{control.observe.choices{len}};
        end
        
    % Populations as projectors onto eigenstates
    case 'prj'
        tise.mat = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            util.disp (['Observable ' int2str(len) ': Populations of eigenstates: ' tise.lab{len}])
            tise.mat{len} = zeros (time.main.n);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                tise.mat{len}(ii,ii) = 1;
            end
        end
        
    % Populations from overlaps with eigenstates
    case 'ovl'
        tise.vec = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            util.disp (['Observable ' int2str(len) ': Overlaps with eigenstates: ' tise.lab{len}])
            tise.vec{len} = zeros (time.main.n,1);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                tise.vec{len}(ii,1) = 1;
            end
        end
        
    otherwise
        util.error('Wrong choice of observable types')
end

%% Save structure 'tise' to disk
util.disp (' ')
util.disp ('Saving all matrix representations in file: tise.mat')
save('tise','tise');

%% Plot matrices

figure(11);
clf;
plot.logo
global plots

% Energy level diagram (diagonal matrix)
subplot(2,2,1)
plot(0:time.main.n-1,real(tise.ham),'o')
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
xlabel('n')
ylabel('E(n)')
title('Energy level diagram')

% System-bath coupling
if isfield (tise, 'sbc')
    subplot(2,2,2)
    surf(tise.sbc)
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    xlabel('n')
    ylabel('m')
    zlabel('\chi(n,m)')
    title('System-bath coupling')
end

% Dipole moments (x)
if isfield (tise, 'd_x')
    subplot(2,2,3)
    surf(tise.d_x)
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    xlabel('n')
    ylabel('m')
    zlabel('\mu_x(n,m)')
    title('Dipole moments (x)')
end

% Dipole moments (y)
if isfield (tise, 'd_y')
    subplot(2,2,4)
    surf(tise.d_y)
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    xlabel('n')
    ylabel('m')
    zlabel('\mu_y(n,m)')
    title('Dipole moments (y)')
end

% Output clock/date/time
util.clock;

end



