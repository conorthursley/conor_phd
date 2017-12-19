% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.


% This script calculates and returns the cross-correlation between two
% calculations.
%
% Both calculations have to be equal in all respects except possibly the grid
% size. The script will load both calculations, convert them to a common grid,
% and calculate the cross correlation. The use-cases are:
%
% 1. Checking for convergence. You run two calculations with different grid
%    parameters. The difference between the absolute value of the cross
%    correlation and 1 gives you a hand-waving size of the numerical error in
%    your calculation.
% 2. Checking that something worked properly. You changed some code or
%    whatever, and want to know precisely, but in a simple way if it works.
%    Just rerun the calculation and compare the wave functions before and
%    after.
%
% As input, give the directories and file names (values of ket.save.dir and
% ket.save.file) for both calculations, and the range of time steps that you
% want to calculate the cross correlation for. If the last parameter is empty,
% we calculate it for all time steps.
%
% The output is the correlation between the two calculations, and the norm of
% both (to exclude small effects).

function [corr, norm1, norm2] = qm_correlation(dir1, file1, dir2, file2, steps)

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

%% Load the calculations, and setup variables
context1 = ket.load(dir1, file1, true);
context2 = ket.load(dir2, file2);

if nargin < 5 || isempty(steps)
	steps = 1:context1.time.main.n;
end

corr = zeros(numel(steps), 1);
norm1 = zeros(numel(steps), 1);
norm2 = zeros(numel(steps), 1);


for index = 1:numel(steps)
	context1 = ket.load(context1, steps(index));
	context2 = ket.load(context2, steps(index));

    wf1 = context1.wf.grid_ND;
    wf2 = context2.wf.grid_ND;

	for m = 1:context1.hamilt.coupling.n_eqs
		norm1(index) = norm1(index) + sum(abs(wf1{m}(:)).^2 .* context1.space.dvr.weight_ND(:));
        norm2(index) = norm2(index) + sum(abs(wf2{m}(:)).^2 .* context2.space.dvr.weight_ND(:));

        [trafo xxx] = ket.transform(wf1{m}, size(wf2{1}), false);
        corr(index)  = corr(index) ...
                + sum(conj(trafo(:)) .* wf2{m}(:) .* context2.space.dvr.weight_ND(:));
    end
end

corr = abs(corr).^2 ./ (norm1 .* norm2);
