% Copyright (C) 2009-2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% Calculate the absorption spectrum of the S1 state only. For this, we have to
% do the following
%
% a) obtain the two lowest torsional states of S0
% b) run a simulation to obtain the eigenstates of the S1 state.
% c) for each eigenstate, calculate the sum of the overlap with the
%    S0 states. This gives the peak height, the eigenergy gives the peak
%    position.
% d) Collect all the peak positions/heights and plot them

qm_setup();
global atomic psi expect space

% a) Obtain the two lowest S0 states. We expect them to lie in 
%    ../1/state_*.dat. If they are not there, throw an error. The grid
%    should also be the same, so we just need to load the values of psi.
%    This is somewhat crude, but should suffice and takes up only few lines
if length(dir('../1/state_*.dat')) ~= 2
    util.error(['The runscript expects the files state_g.dat and state_g.dat'...
                ' in directory ../1 relative to current working path!']);
end

% b) and c). This has to be done twice because for some reason, the
% diagonalisation LAPACK routine is not aware of the symmetry of our
% problem.
for symmetry = ['g' 'u']
    % First, run the simulation with the given symmetry
	qm_init(symmetry);
    qm_bound();

    % an offset to store the ungerade states after the gerade ones.
    if symmetry == 'u'
        offset = psi.eigen.stop - psi.eigen.start + 1;
    else
        offset = 0;
    end

    % load the corresponding state
    mystate = load(strcat(['../1/state_' symmetry '.dat']), '-ascii');
        
    % Finally the heart of the routine: Calculate the spectrum; the dipole moment
    % is constant 1 a.u.
    for step = psi.eigen.start+1:psi.eigen.stop+1
        ket.eigen(step);

        positions(step - psi.eigen.start + offset) = expect.tot.tot(step) * atomic.w.cm_1;

        heights(step - psi.eigen.start + offset)  = ...
                abs( sum(conj(psi.dvr.grid_ND{1}(:)) .* mystate .* space.dvr.weight_ND(:)) ).^2;
    end
end

qm_cleanup();

% Smooth out the peaks by a Gaussian with 3 cm^-1 FWHM, and migrate the data
% on an equally spaced grid.
energies = 25700:0.1:27000;
peaks    = zeros(size(energies));

for ii = 1:length(energies)
    factors = exp(-(positions-energies(ii)).^2 / (2 * 3/sqrt(8 * log(2))).^2);
    peaks(ii) = sum(heights(:) .* factors(:));
end

% Plot the resulting spectrum (should it be smoothed in-between?)
figure(1)
clf
plot(energies, peaks, 'LineWidth', 2);
xlabel('E (cm^-1)')
saveas(gcf, 'excite.jpg')
