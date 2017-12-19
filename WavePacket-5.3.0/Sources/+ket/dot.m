% Calculate inner product <bra|ket> by DVR (quadrature) methods

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.

function retval = dot (bra, ket)
global space hamilt

retval = 0;

for m = 1:hamilt.coupling.n_eqs
    retval = retval + sum(conj(bra{m}(:)) .* ket{m}(:) .* space.dvr.weight_ND(:));
end

end

