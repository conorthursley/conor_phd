%------------------------------------------------------------------------------
%
% For given indices of first (time.main.start) and last (time.main.stop) time step
% and for given size of timesteps (time.main.delta) and given number of substeps
% (time.sub.n) per time step, this function calculates the total simulation
% time  (time.main.total) and size of sub steps (time.sub.delta) 
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2010-2011 Ulf Lorenz
%
% see the README file for license details.

function timesteps 
global time hamilt

%% Main steps: by default starting from zero
if ~isfield (time.main,'start')
    time.main.start = 0;
end
time.main.n      = abs(time.main.stop - time.main.start) + 1;              % number of time steps
time.main.total  = time.main.delta * ( time.main.n - 1 );                  % total simulation time
time.main.grid   = (time.main.start + [ 0 : time.main.n-1 ]') * time.main.delta; % Column(!) vector


%% Substeps: default number is one
if ~isfield (time,'sub') || ~isfield (time.sub,'n')
    time.sub.n = 1;
end
time.sub.delta   = time.main.delta / time.sub.n;        % size of substeps
time.sub.grid    = time.main.start*time.main.delta ...
        + [0:(time.main.n-1)*time.sub.n]' * time.sub.delta;

util.disp ('*******************************************************')
util.disp ('Temporal discretization: Constant steps and sub-steps')
util.disp ('*******************************************************')
util.disp ( [ 'Index of first step  : ' int2str(time.main.start) ] )
util.disp ( [ 'Index of last step   : ' int2str(time.main.stop)  ] )
util.disp ( [ 'Number of steps      : ' int2str(time.main.n-1)  ] )
util.disp ( ' ' )
util.disp ( [ 'Size of steps        : ' num2str(time.main.delta)  ] )
util.disp ( [ 'Total time           : ' num2str(time.main.total)   ] )
util.disp ( ' ' )
util.disp ( [ 'Number of sub-steps  : ' int2str(time.sub.n)  ] )
util.disp ( [ 'Size of sub-steps    : ' num2str(time.sub.delta)  ] )
util.disp ( ' ' )

% hamilt.range.delta will be only set/used when this function
% is called from qm_propa but not from qm_control, qm_optimal
if isfield (hamilt,'range')
    if isfield (hamilt.range,'delta')
        time.main.alpha  = time.main.delta  * hamilt.range.delta / 2;
        time.sub.alpha   = time.sub.delta   * hamilt.range.delta / 2;
        util.disp ( [ 'Spectral range of Hamiltonian : ' num2str(hamilt.range.delta)  ] )
        util.disp ( [ 'Kosloff alpha (main steps)    : ' num2str(time.main.alpha)  ] )
        util.disp ( [ 'Kosloff alpha (sub steps)     : ' num2str(time.sub.alpha)  ] )
        util.disp ( ' ' )
    end
end
