% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


% Run the relaxation, and copy the resulting relaxed wavefunction to the
% other run directories so that it can be used as input for the interpolating
% initial wave function
qm_setup();
qm_init();
qm_propa();

global space psi

output = cat(2, space.dvr.grid_ND{1}, psi.dvr.grid_ND{1});
dlmwrite('../2/wav_1.dat', output, ' ');

qm_cleanup();
