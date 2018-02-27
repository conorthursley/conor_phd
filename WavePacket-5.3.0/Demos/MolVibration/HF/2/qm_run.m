% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% Calculate and plot the pulse shape and frequency, and the population of
% several eigenstates of a Morse oscillator.
%
% This script aims at reproducing figs 1, 2 of Phys.Rev.Lett 65:2355. For this,
% we first run a calculation that performs the actual propagation in question,
% then we misuse various features of WavePacket to obtain what we want.

%% 1. Run the calculation. It should save the results in the current directory,
%     in a file 'HF'.
qm_setup();
qm_init();
qm_propa();
qm_cleanup();

%% 2. Get a handle for loading the results of the calculation
context = ket.load('.', 'HF');
init.info('rundemo');


%% 3. Initialise the eigenstates of the Morse potential
% The function wav_morse already calculates eigenstates of the Morse potential.
% Though it is meant to do so for the purpose of setting up the initial wave
% function, we can easily tweak it to give us the Morse eigenfunctions.
%
% Prior knowledge: The Morse potential of HF supports 24 bound states. The
% calculation has been set up with the Morse ground state; we only need to
% adjust the quantum number.

global psi

morse = cell(24, 1);
for ii = 1:24
    psi.dof{1}.n_q = ii-1;

    morse{ii} = wav.morse(1);
end


%% 4. Obtain the shape of the electric field, and the frequency. We calculate
%     them by some code snipplets copied from util.efield.m. Note that the
%     frequency calculation only works like this because the pulses are not
%     expected to overlap, otherwise the frequencies of different pulses would
%     add up to some unphysical quantity.

global time

% First, the envelope
total_envelope = zeros(size(time.main.grid));
for ii = 1:time.efield.n_pulse
    tt0 = time.main.grid - time.efield.delay(ii);

    switch lower(time.efield.shape(ii,:))
        case 'recta'
            envelope = ones(size(tt0));
        case 'inter'
            envelope = interp1(time.efield.times{ii}, time.efield.fields{ii}, ...
                    tt0, time.efield.method, 0);
    end

    where = find( abs(tt0) > time.efield.length(ii)/2 );
    envelope(where) = 0;
    total_envelope = total_envelope + envelope;
end

% Then the frequency.
all_frequencies = zeros(size(time.main.grid));
for ii = 1:time.efield.n_pulse
    tt0 = time.main.grid - time.efield.delay(ii);

    frequency = time.efield.frequ(ii) * ones(size(tt0));
    frequency = frequency + time.efield.linear(ii) * tt0 ...
                + time.efield.quadratic(ii)/2 * tt0.^2;

    % somewhat annoying again: we wish to define the frequency as w(t)
    % in cos ( w(t) * t ). However, we have defined it here as w(tt0) as
    % in cos( w(tt0) * tt0 + phi ). In general, these two definitions give
    % different results, so we have to transform. Equalling both definitions
    % gives for t ~= 0:
    % w(t) = [ w(tt0) * tt0 + phi ] / t
    frequency(2:end) = ( frequency(2:end) .* tt0(2:end) + time.efield.phase(ii) ) ...
            ./ time.main.grid(2:end);
    frequency(1) = frequency(2);        % prior knowledge

    where = find( abs(tt0) > time.efield.length(ii)/2 );
    if strcmp(lower(time.efield.shape(ii,:)), 'inter')
        support = time.efield.times{ii};
        where = find( tt0 > max(support(:)) | tt0 < min(support(:)));
    end
        
    frequency(where) = 0;
    all_frequencies = all_frequencies + frequency;
end

% normalize to the transition frequency between the lowest two levels of the
% Morse potential. These happen to be correlated with the main time steps,
% which makes this fairly simple.
all_frequencies = all_frequencies * time.main.delta / (2*pi);


%% 5. Calculate the time-dependent projections on the Morse eigenstates.
pdiss = ones(time.main.n, 1);
projection = cell(size(morse));
for ii = 1:length(projection)
    projection{ii} = zeros(time.main.n, 1);
end

for timestep = 1:time.main.n
    context = ket.load(context, timestep);

    for ii = 1:length(projection)
        projection{ii}(timestep) = projection{ii}(timestep) ...
                + abs( sum(conj(context.wf.grid_ND{1}(:)) .* morse{ii}(:) .* context.space.dvr.weight_ND(:)) )^2;
    end
end

for ii = 1:length(projection)
    pdiss = pdiss - projection{ii};
end

%% 6. Now visualize everything.
timesteps = time.main.start:time.main.stop;

% The projections
figure(1);
clf;
hold on;

plot(timesteps, projection{1}, 'k--', ...
    timesteps, projection{2}, 'k-', ...
    timesteps, projection{9}, 'k-', ...
    timesteps, projection{15}, 'k--', ...
    timesteps, pdiss, 'k-', ...
    'LineWidth', 2);
xlabel('time (cycles)');
ylabel('probabilities');
axis( [0 500 0 1.1] );

hold off;
saveas(gcf, 'projections.jpg');

% The envelope and frequency
figure(2);
clf;
hold on;

plot(timesteps, total_envelope, 'k-', ...
    timesteps, all_frequencies, 'k--', ...
    'LineWidth', 2);
xlabel('time (cycles)');
axis([0 500 0 1.2]);

hold off;
saveas(gcf, 'pulse.jpg');

qm_cleanup();
