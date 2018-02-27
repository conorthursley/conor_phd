% Copyright (C) 2008-2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


%---------------------------------------------------------------
% This script intends to reproduce Figure 3
% of [ J.Chem.Phys 110:3870].
%
% Note that the figure in the paper contains some errors:
% * the average over all M should start at <cos^2> = 0.33.
% * the curve labelled M = 4 is probably M=5.
%---------------------------------------------------------------

qm_setup();

%% Do the propagation for all 6 values of M and
%% save the resulting expectation values.
global expect time

average = cell(6,1);

for mval = 0:5
	qm_init(mval);
	qm_propa();
	average{mval+1} = expect.tot.amo;
end

averagetot = average{1}/11 + (average{2} + average{3} + average{4} + average{5} + average{6}) *2/11;

%% Plot it in the same way as in the paper.
figure(1);
clf;
plot(time.main.grid, averagetot, 'k-', ...
     time.main.grid, average{1}, 'k--', ...
     time.main.grid, average{3}, 'k-.', ...
     time.main.grid, average{5}, 'k:', ...
     'LineWidth', 2);
xlabel('t');
ylabel('<cos^2 \theta>');
saveas(gcf, 'alignment.jpg');
