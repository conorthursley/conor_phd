% Copyright (C) 2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

% TODO: Does not work right now, needs an update.


% This is a demo script for doing optimal control theory with WavePacket. We use
% the optimal control scheme of Kosloff, Rice, Tannor et al. to be found in 
% Chem. Phys. 139:201 (1989), and apply this to the problem of optimal HOD bond
% cleavage.
%
% Due to the way WavePacket is written, and because of the goal of making this
% script easy to parse, there are a few things that are not optimal, and which I
% would point out here.
% 
% 1. The execution speed could be easily improved by about 1/3.
%
%    Kosloff et al. start with a trial electric field, propagate their wave
%    packet forward in time, and then back-propagate the final wave function and
%    the desired target state. Obviously, if you had saved the wave function
%    during the forward-propagation, propagating it backwards in time again
%    would just be silly, since you get nothing new.
%
%    The reason this is still done here is that it saves some ugly index
%    manipulation when you later load the propagated wave function and the
%    result of the back-propagation of the target state. If you understand what
%    this text means and are confident that this is not difficult to overcome
%    (it is not), feel free to change a few lines of code and get a major speed
%    boost.
%
% 2. The script is not as straight-forward and fast as it could be.
%
%    WavePacket as of now only supports the propagation of _one_ wave function
%    at a time. After each main time step, expectation values are calculated
%    etc., and the wave funtion is saved. Since we need, however, two wave
%    functions (the propagated state and the propagated target states), we can
%    overcome this limitation by running two propagations, saving the results,
%    and then loading everything again to calculate the overlaps.
%
%    This procedure sounds not too complicated, but leads to a serious amount of
%    overhead in the script. Also, for example, we need to save the wave
%    function at a relatively fine-grained time stepping, which leads to some
%    large consumption of disk space and makes the whole thing a bit slow.
%    However, these shortcomings can only be fixed by a major rewrite of the
%    code base (e.g., adapting the qm_propa script), which is outside the scope
%    of this demo.

% params is used to supply some parameters to the input script
% options holds all global settings for the script
global atomic params options;

global psi space


%% Parameter input
%
% Set all sorts of options for the script here.


% 1. Define the initial state.
%
%    This setting is copied to psi.init whenever we propagate our initial state
%    with one of the trial electric fields.
options.init.dof{1}.handle = @wav.gauss;
options.init.dof{1}.width  = 1/sqrt(4*23.18);
options.init.dof{1}.pos_0  = 1.806;
options.init.dof{1}.mom_0  = 0;

options.init.dof{2}.handle = @wav.gauss;
options.init.dof{2}.width  = 1/sqrt(4*18.84);
options.init.dof{2}.pos_0  = 1.833;
options.init.dof{2}.mom_0  = 0;

options.init.coeffs        = [1 0];
options.init.representation= 'dia';

% 2. Define the projection operator.
options.amo.handle = @projector;

% 3. Define the time steps (fast for the forward propagation and slow for the 
%    calculation of the electric field) and some laser parameters.
options.time.fast.main.start = 0;
options.time.fast.main.stop  = 10;
options.time.fast.main.delta = 225;
options.time.fast.sub.n      = 225;

% Note the negative time, because we always use units of time.main.delta (which
% is negative), i.e., time.main.start * time.main.delta should be the final
% time of the forward propagation.
options.time.slow.main.start = -2250;
options.time.slow.main.stop  = 0;
options.time.slow.main.delta = -1;
options.time.slow.sub.n      = 1;

% 4. The initial guess of the electric field, and the setup of the electric field
%    We consider the 56 TW/cm^2 peak electric field case.
options.efield.init.dressed = false;
options.efield.init.shape   = 'gauss';
options.efield.init.ampli   = sqrt(56/atomic.I.TW_cm2);
options.efield.init.delay   = 400;
options.efield.init.frequ   = 0.2738;
options.efield.init.phase   = 0;
options.efield.init.polar   = 0;
options.efield.init.fwhm    = 5/atomic.t.fs;

options.efield.load.dressed = false;
options.efield.load.shape   = 'inter';
options.efield.load.complex = true;
options.efield.load.delay   = 0;
options.efield.load.frequ   = 0;
options.efield.load.phase   = 0;
options.efield.load.polar   = 0;
options.efield.load.method  = 'spline';
% the file name is determined on the fly.

% keep the total energy constant
options.efield.energy = 1/2 * options.efield.init.ampli^2 * options.efield.init.fwhm ...
                        * sqrt(pi) / (2 * log(2));

% 4. Various file names and other parameters
options.save.dir    = pwd;
options.save.forward= 'forward';
options.save.init   = 'initial';
options.save.amo   = 'projected';

options.efield.file = 'field';     % full filename is "${file}_${iteration}.dat"

options.iterations  = 15;


%% Now run the calculation

for iter = 1:options.iterations
    % this prints out the text in the red signal color, so it is visible among
    % all the output from the propa scripts.
    fprintf(2, 'Iteration %i', iter);

    % A) Propagate the initial state forward in time, and save it. This is done
    %    by a simple call to qm_propa; what takes code is setting up the
    %    parameters.
    params.time = options.time.fast;
    if iter == 1
        params.time.efield = options.efield.init;
    else
        params.time.efield = options.efield.load;
        params.time.efield.file = {[options.efield.file '_' num2str(iter-1) '.dat']};
    end

    params.psi.init        = options.init;
    params.psi.save.export = true;
    params.psi.save.dir    = options.save.dir;
    params.psi.save.file   = options.save.forward;

    params.plot = false;

    qm_propa;


    % B) Propagate the initial state backwards again
    params = [];

    params.time = options.time.slow;
    if iter == 1
        params.time.efield = options.efield.init;
    else
        params.time.efield = options.efield.load;
        params.time.efield.file = {[options.efield.file '_' num2str(iter-1) '.dat']};
    end

    params.psi.corr.handle = @wav.load;
    params.psi.corr.dir    = options.save.dir;
    params.psi.corr.file   = options.save.forward;
    params.psi.init.representation = 'dia';
    params.psi.save.export      = true;
    params.psi.save.dir         = options.save.dir;
    params.psi.save.file        = options.save.init;

    params.plot = false;

    qm_propa;

    % C) Project the final wave function on the desired exit channel, and
    %    propagate the result backwards in time. Now this requires some manual
    %    intervention.

    %  We first manually assemble a custom save file from which we later load
    %  the "initial state".
    context = ket.load(options.save.dir, options.save.forward, true);
    context = ket.load(context, 2);
    
    init.info('tmp');
    psi.save.dir = options.save.dir;
    psi.save.file= 'tmp';

    % fill up the save file with dummy wave functions.
    for tt = 1:context.time.main.n-1
        ket.save(tt);
    end

    % the last wave function is the only one that is used by wav_load() anyway.
    psi.dvr.grid_ND = feval(options.amo.handle, context.wf.grid_ND);
    ket.save(context.time.main.n);
    fclose('all');


    params = [];

    params.time = options.time.slow;
    if iter == 1
        params.time.efield = options.efield.init;
    else
        params.time.efield = options.efield.load;
        params.time.efield.file = {[options.efield.file '_' num2str(iter-1) '.dat']};
    end

    params.psi.corr.handle = @wav.load;
    params.psi.corr.dir    = options.save.dir;
    params.psi.corr.file   = 'tmp';
    params.psi.init.representation = 'dia';
    params.psi.save.export      = true;
    params.psi.save.dir         = options.save.dir;
    params.psi.save.file        = options.save.amo;

    params.plot = false;

    qm_propa;


    % D) Now we can go on and calculate the overlap between the functions, on
    %    which we base the new electric field.
    context_init = ket.load(options.save.dir, options.save.init);
    context_amo  = ket.load(options.save.dir, options.save.amo);

    efield = zeros(context_init.time.main.n, 2);
    efield(:, 1) = context_init.time.main.grid;

    for step = 1:size(efield, 1)
        context_init = ket.load(context_init, step);
        context_amo  = ket.load(context_amo, step);

        wf_init = context_init.wf.grid_ND;
        wf_amo  = context_amo.wf.grid_ND;
        trans_dipole = context_init.hamilt.d_x.grid_ND{1,2};

        efield(step, 2) = 1i * (...
                sum( conj(wf_amo{2}(:)) .* trans_dipole(:) .* wf_init{1}(:) .* space.dvr.weight_ND(:)) ...
                - sum( conj(wf_init{2}(:)) .* trans_dipole(:) .* wf_amo{1}(:) .* space.dvr.weight_ND(:)) );
    end


    % E) Normalize the electric field and write it out to a file
    efield(:, 2) = efield(:, 2) * sqrt(options.efield.energy) ...
             / sqrt( sum(abs(efield(:,2)).^2) * abs(context_init.time.main.delta) );

    filename = [options.efield.file '_' num2str(iter) '.dat'];
    dlmwrite(filename, cat(2, efield(:, 1), real(efield(:, 2)), imag(efield(:, 2))), ...
            'precision', 10, 'delimiter', ' ');
end
