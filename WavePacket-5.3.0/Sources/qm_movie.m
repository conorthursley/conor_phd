%--------------------------------------------------------------------------
%
% Creates movie output from previous qm_bound or qm_propa calculations
%
% qm_movie(savedir, savefile) loads the saved calculation in the file
% "savefile" found in directory "savedir" and plots all graphs.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2010-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_movie(savedir, savefile)

global info plots time;

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Provide default values for missing input arguments
if nargin<2
    savefile = 'WavePacketSave';
end

if nargin<1
    savedir = pwd;
end

util.disp('-------------------------------------------------');
util.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_movie')
util.disp('Make movie from saved calculation in file        ');
util.disp(['     ' savefile]);
util.disp('residing in directory                            ');
util.disp(['     ' savedir]);
util.disp('-------------------------------------------------');
util.disp(' ');

% Load the general data of the saved calculation.
% Returns a context that we can later use to 
% load the wave function at specific time steps.
% Sets most global variables, except for: info, plots
context = ket.load(savedir, savefile, true);

% the standard plots are enabled by default
plots.density.on = true;
plots.expect.on = true;
plots.density.export.on = true;
plots.expect.export.on = true;

% copy the program name from the saved file; this is needed in some checks
info.program = context.info.program;

% reset expectation values (calculating them live should not waste so much
% computing time, and it gives a visible progress bar)
init.expect();

% Main loop: For each time step, load the wave function from the calculation
for step  = 1:time.main.n
	context = ket.load(context, step, true);

	ket.adiabatic('dia2adi')

	ket.expect(step);

	util.logging(step)

	plot.wigner(step);
	plot.psi(step);

	ket.adiabatic('adi2dia')
end

% Output clock/date/time
util.clock;

end
