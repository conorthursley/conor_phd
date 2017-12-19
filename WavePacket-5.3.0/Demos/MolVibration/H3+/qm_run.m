% Calculates the lowest three vibrational states of H3+. See work by
% W. Meyer, P. Botschwina, P. Burton, J. Chem. Phys. 84:891 (1986)
% According to Table VIII of that reference (doi:10.1063/1.450534),
% we should see the eigenstates at 4345, 6861, 7530 cm^-1.
% Due to convergence problems, our results may deviate by about
% 10-20 cm^-1.
%
% This script uses propagation in imaginary time to get the
% eigenstates. It will first calculate the groundstate, then the 
% first excited state, then the second one. Each of the latter two
% calculations uses a hacked version of the ket.cheby_imag propagator
% that propagates only in the space perpendicular to the lower 
% eigenstates. All three propagation (subdirs gs, exc1, exc2) must
% have the same grid parameters!

global atomic gs energy

%% First, do the calculation of the groundstate and save it
cd gs
qm_setup(); qm_init(); qm_propa(); qm_cleanup();

global psi space
gs = psi.dvr.grid_ND{1};

% Calculate and store the energy
energy = zeros(3,1);

ket.hamilt(0,0,0);
energy(1) = sum( conj(psi.dvr.new_ND{1}(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );


%% Calculate the first excited state and store it as well
cd ../exc1
qm_setup(); qm_init(); qm_propa(); qm_cleanup();

global psi space
temp = gs;
gs = cell(2,1);
gs{1} = temp;
gs{2} = psi.dvr.grid_ND{1};

% Calculate and store the energy
ket.hamilt(0,0,0);
energy(2) = sum( conj(psi.dvr.new_ND{1}(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );


%% Calculate the second excited state
cd ../exc2
qm_setup(); qm_init(); qm_propa(); qm_cleanup();

% Calculate and store the energy
global psi space

ket.hamilt(0,0,0);
energy(3) = sum( conj(psi.dvr.new_ND{1}(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );


%% Output
cd ..
energy = abs(energy) * atomic.w.cm_1;    % conversion to cm^-1

disp( 'Energies calculated (cm^-1):' )
disp( ['0: ' num2str(energy(1))] )
disp( ['1: ' num2str(energy(2))] )
disp( ['2: ' num2str(energy(3))] )
