% Calculate matrix elements ("sandwiches") of an operator <bra|operator|ket> 
% by means of DVR (quadrature) methods


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%
% see the README file for license details.

function entry = sandwich (bra, operator, ket)
    global space hamilt

    entry = 0;

    for m = 1:hamilt.coupling.n_eqs
        % diagonal couplings count only once
        if ~isempty(operator{m,m})
            entry = entry + sum(conj(bra{m}(:)) .* operator{m,m} .* ket{m}(:) .* space.dvr.weight_ND(:));
        end

        % the offdiagonal coupling maps both off-diagonal elements on one term only
        for n = m+1:hamilt.coupling.n_eqs
            if ~isempty(operator{m,n})
                entry = entry ...
                        + sum(conj(bra{m}(:)) .* operator{m,n}(:) .* ket{n}(:) .* space.dvr.weight_ND(:)) ...
                        + sum(conj(bra{n}(:)) .* conj(operator{m,n}(:)) .* ket{m}(:) .* space.dvr.weight_ND(:));
            end
        end
    end
    
end

