% Copyright (C) 2008-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


%---------------------------------------------------------------
% This script intends to reproduce the first graph in Figure 1
% of [J.Chem.Phys 110:3870].
%---------------------------------------------------------------

% Clean workspace and set up the basic propagation.
qm_setup();

%% Do the propagation for all three values of Delta omega and
%% save the resulting expectation values.
global time expect

qm_init(100);
qm_propa();
expect1 = expect.tot.amo;
mom1    = expect.ind.fbr2{1};

qm_init(400);
qm_propa();
expect2 = expect.tot.amo;
mom2    = expect.ind.fbr2{1};

qm_init(900);
qm_propa();
expect3 = expect.tot.amo;
mom3    = expect.ind.fbr2{1};


%% Plot it in the same way as in the paper.
figure(1);
clf;
plot(time.main.grid, expect1, 'k-', ...
     time.main.grid, expect2, 'k--', ...
     time.main.grid, expect3, 'k:', ...
     'LineWidth', 2);
xlabel('t');
ylabel('<cos^2 \theta>');
saveas(gcf, 'alignment.jpg');

figure(2);
clf;
plot(time.main.grid, mom1, 'k-', ...
     time.main.grid, mom2, 'k--', ...
     time.main.grid, mom3, 'k:', ...
     'LineWidth', 2);
xlabel('t');
ylabel('<J^2>');
saveas(gcf, 'j2.jpg');

qm_cleanup();
