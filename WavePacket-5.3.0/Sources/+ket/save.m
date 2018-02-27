% ------------------------------------------------------------------------
%
% Stores the wavefunction in an internal variable and finally saves it
% in a MATLAB® formatted binary file (MAT-file, file name extension .mat) 
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%
% see the README file for license details.

function save ( step )

global expect hamilt info plots psi space time uncert

% index is the index in save_wf where we put the next WF.
% save_wf is a cell array with each element containing a "time" and a "grid_ND"
%    information
persistent store index save_wf


%% Initialize everything
if step == 1
    % check if the wavefunction shall be stored at all.
    util.disp('****************************************');
    if ~isfield(psi, 'save') || ~isfield(psi.save, 'export') || ~psi.save.export
        util.disp('Not saving the wave function')
        store = false;
    else
        store = true;
        index = 1;
        
        if ~isfield (psi.save,'dir') % default directory: where you are now
            psi.save.dir = pwd;
        end
        
        if ~isfield (psi.save,'file') % default file name template
            psi.save.file = 'WavePacketSave';
        end
        
        util.disp(['Wave function saved to dir "' psi.save.dir '",file(s) "' psi.save.file '"'])

        % Determine the stepsize; By default, we write data in 500 MB chunks,
        % this limitation is there to avoid blasting the system memory.
        if ~isfield(psi.save, 'mem')
            psi.save.mem = 500 * 2^20;
        end
        mem_per_step = numel(space.dvr.grid_ND{1}) * 16 * hamilt.coupling.n_eqs;
        psi.save.stepsize = min(int32(ceil(psi.save.mem / mem_per_step)), ...
                time.main.n+1);

        if isfield(psi.save, 'step') && psi.save.stepsize > psi.save.step
            psi.save.stepsize = psi.save.step;
        end

        % Create a cell array for the saved wavefunction
        save_wf = cell(psi.save.stepsize, 1);
    end
    util.disp('****************************************');
    util.disp(' ');
    util.disp(' ');
end


%% If we do not store the wavefunction, exit immediately
if ~store
    return;
end


%% Store the wavefunction
save_wf{index} = psi.dvr.grid_ND;
index = index + 1;


%% Write out the WF and possible additional data.
if index > psi.save.stepsize || step == time.main.n
	if ~isdir(psi.save.dir)
		mkdir(psi.save.dir);
	end

	file_index = int32(floor( (step-1)/double(psi.save.stepsize) ));
    filename = fullfile(psi.save.dir, strcat(psi.save.file, '_', int2str(file_index), '.mat'));
    save(filename, 'save_wf', '-mat');

    index = 1;
    save_wf = cell(psi.save.stepsize, 1);
end

if step == time.main.n
    filename = fullfile(psi.save.dir, strcat(psi.save.file, '.mat'));
    save(filename, 'expect', 'hamilt', 'info', 'plots', 'psi', ...
                   'space', 'time', 'uncert', '-mat');
end
