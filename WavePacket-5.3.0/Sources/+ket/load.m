% ------------------------------------------------------------------------
%
% Loads the wavefunction and the setup previously stored with ket.save
% from an external file.
%
% This function takes two variables as input, and returns a single variable as
% output. There are two different ways to use this function.
%
% 1. If the input consists of two strings, the first parameter is interpreted as the
%    path, the second as the filename of the saved calculation to load. The
%    format is the same as the one used for saving the wavefunction. The function
%    then returns a structure that contains all the global variables from the
%    run.
%
%    If you call
%
%    context = ket.load('mydir', 'myfile');
%
%    you can then access context.space, context.hamilt etc. to get the
%    information that the run was set up with.
%
%    If a third parameter is supplied (boolean with true value), this call
%    also overwrites all global variables with those from the saved calculation.
%    One global variable that is never set is info.
%
%    The context structure contains two additional fields. The first one is
%    context.cache, which holds the wave functions from one complete save file
%    for more efficient access. The other is context.wf, where the wave
%    functions are stored (see below).
%
% 2. If you supply such a handle context, and an integer t, the function loads the
%    wavefunction from the file specified within context, at timestep t
%    (starting at 1 for the initial time), and sets context.wf.grid_ND to give
%    the wave function as a grid. If t is an array of integers, context.wf is a
%    cell array whose n-th element is the wavefunction at the timestep t(n).
%
%    As an example, if you call
%
%    context = ket.load(context, [1 2]);
%
%    you can calculate the norms of the first state at the first two timesteps via
%
%    norm_at_1 = sum(abs(context.wf{1}.grid_ND{1}(:)).^2 .* context.space.weight_ND(:));
%    norm_at_2 = sum(abs(context.wf{2}.grid_ND{1}(:)).^2 .* gv.space.weight_ND(:));
%
%    If a third parameter is supplied with boolean value of true, this call also
%    stores the wavefunction (last wave function if multiple timesteps are
%    given) in psi.dvr.grid_ND.
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%
% see the README file for license details.

function context = load(arg1, arg2, glob)
global expect hamilt psi space uncert time;

%% Two strings as parameters => set up the loading stuff
if isa(arg1, 'char') && isa(arg2, 'char')
    filename = fullfile(arg1, strcat(arg2, '.mat'));
    context   = load(filename, '-mat');
    context.psi.load.dir  = arg1;
    context.psi.load.file = arg2;
	context.cache = [];
	context.wf    = [];
    % if the additional parameter glob is true
    % set most global variables, except for: info, plots
    if nargin == 3 && glob == true
        expect = context.expect;
        hamilt = context.hamilt;
        psi    = context.psi;
        space  = context.space;
        uncert = context.uncert;
        time   = context.time;
    end
%% Otherwise, first parameter is the context, second is the time steps
elseif isfield(arg1, 'psi') && isa(arg2, 'numeric')
    context   = load_from_context(arg1, arg2);
    % if the additional parameter glob is true, set psi.dvr.grid_ND as well
    if nargin == 3 && glob == true
		if numel(arg2) > 1
	        psi.dvr.grid_ND = context.wf{numel(arg2)}.grid_ND;
		else
			psi.dvr.grid_ND = context.wf.grid_ND;
		end
    end
%% Everything else is an error
else
    util.error(['Tried to load with incorrect data types. Arguments must either both ' ...
        'be strings, or a context and an integer array.']);
end

end % psi.load


function context = load_from_context(context, times)
    % some checks
    if any(times > context.time.main.n) || any(times < 1)
        util.error('time index out of range')
    end

    if numel(times) > 1
        context.wf = cell(numel(times), 1);
    end

    for ii = 1:numel(times)
        step = times(ii);
        
        % Calculate indices for the file suffix and the position in the file.
        ind_file = int32(floor( (step-1)/double(context.psi.save.stepsize) ));
        ind = step - ind_file*context.psi.save.stepsize;

        % Check if we have cached this file, otherwise load and cache it
        if isempty(context.cache) || ind_file ~= context.cache.file_index
            filename = fullfile(context.psi.load.dir, ...
                    strcat(context.psi.load.file, '_', int2str(ind_file), '.mat'));
            mydata = load(filename, '-mat');

            context.cache.file_index = ind_file;
            context.cache.data = mydata.save_wf;
        end

        % Assign the return value
        if numel(times) == 1
            context.wf.grid_ND = context.cache.data{ind};
        else
            context.wf{ii}.grid_ND = context.cache.data{ind};
        end
    end
end
