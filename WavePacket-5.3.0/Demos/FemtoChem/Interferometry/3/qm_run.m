% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% Calculate the pump-probe signal for two-photon interferometry.
% We vary the timesteps between pump and probe in intervals of 10 fs,
% For each delay time, we run a calculation that loads the wavefunction
% from the pump-only simulation, and propagates it for 200 fs. Then, we
% calculate the observable. In our case, this is just the population of
% the 2^1\Pi_g state after the probe.
%
% Note that the reference (J.Cem.Phys. 100:5448) discusses in length the
% population caused by different processes (two photons taken from the
% pump pulse, probe pulse, or one from each). To somewhat disentangle the
% contributions to the final population, we propagate the pump result for
% the Sigma_g, Sigma_u, and Pi_g state separately. For non-overlapping
% pulses, this gives approximately the three wavefunctions psi.{11},
% psi.{22}, and psi.{21}. We can then overlap them as we want, and obtain
% the single contributions. See the reference for more details.
%
% The time delay is passed to the qm_init.m script via the global variable
% delay, the initial coefficients via init.coeffs.
pops.tp1    = zeros(201,1);     % ~2-photon absorbtion from pump pulse
pops.tp2    = zeros(201,1);     % ~2-photon absorbtion from probe pulse
pops.pp     = zeros(201,1);     % pump-probe signal
pops.tpi    = zeros(201,1);     % two-photon interference
pops.pptpi  = zeros(201,1);     % the pump-probe two-photon interference
pops.tpppi  = zeros(201,1);     % the two-photon pump-probe interference
ii = 1;

% for practical reasons, we use the delay time in units of 10 fs.
qm_setup();
global psi space

for delay = 0:1:200
    % get the three wavefunctions
    qm_init(delay, 1);
    qm_propa();
    psi.d22 = psi.dvr.grid_ND{3};

    qm_init(delay, 3);
    qm_propa();
    psi.d11 = psi.dvr.grid_ND{3};

	qm_init(delay, 2);
	qm_propa();
    psi.d21 = psi.dvr.grid_ND{3};

    % calculate all the single terms
    pops.tp1(ii) = sum(abs(psi.d11(:)).^2 .* space.dvr.weight_ND(:));
    pops.tp2(ii) = sum(abs(psi.d22(:)).^2 .* space.dvr.weight_ND(:));
    pops.pp(ii)  = sum(abs(psi.d21(:)).^2 .* space.dvr.weight_ND(:));

    pops.tpi(ii)   = 2 * sum( real(conj(psi.d11(:)) .* psi.d22(:)) .* space.dvr.weight_ND(:));
    pops.pptpi(ii) = 2 * sum( real(conj(psi.d21(:)) .* psi.d22(:)) .* space.dvr.weight_ND(:));
    pops.tpppi(ii) = 2 * sum( real(conj(psi.d21(:)) .* psi.d11(:)) .* space.dvr.weight_ND(:));

    ii = ii + 1;
end

qm_cleanup();

pops.total = pops.tp1 + pops.tp2 + pops.pp + pops.tpi + pops.pptpi + pops.tpppi;

figure(1);
clf;
plot(0:10:2000, pops.total, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('total Population of 2^1\Pi_g');
saveas(gcf, 'total.jpg');

clf;
plot(0:10:2000, pops.tp1 + pops.tp2, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (both photons from one pulse)');
saveas(gcf, 'two-photon.jpg');

clf;
plot(0:10:2000, pops.pp, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (Pump-Probe signal)');
saveas(gcf, 'pp.jpg');

clf;
plot(0:10:2000, pops.tpi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (two-photon interference)');
saveas(gcf, 'tpi.jpg');

clf;
plot(0:10:2000, pops.pptpi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (interference between pump-probe and probe-probe)');
saveas(gcf, 'pptpi.jpg');

clf;
plot(0:10:2000, pops.tpppi, 'LineWidth', 2);
xlabel('delay time (fs)');
ylabel('Population of 2^1\Pi_g (interference between pump-pump and pump-probe)');
saveas(gcf, 'tpppi.jpg');
