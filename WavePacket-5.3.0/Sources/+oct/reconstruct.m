% Reconstructs the raw density operator (in matrix form) in full dimensionality
% to be used for correlation measures, i.e., undo balancing (and truncation)
%
% The input is the current column vector *x*, the equation of motion *eom*
% ('tdse' or 'lvne') determines whether state *x* represents a wave function 
% or a density operator), and the balancing/truncation *scheme*.

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2012 Burkhard Schmidt, Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

function rho = reconstruct(x, eom, scheme)

global control matrix;

% Undo balancing (and truncation)
if scheme > 0
    x = matrix.T * x;
end

switch lower(eom)
    case 'lvne'
        
        % Dimensionality
        dim = round(sqrt(numel(x)));
        
        % undo diagonal first order
        if strcmp('df', control.lvne.order)
            U = util.cw2df(dim);
            x = U.' * x;
        end

        rho = reshape(x, dim, dim);

    case 'tdse'

        % Undo balancing/truncation
        if scheme > 0
            x = matrix.T * x;
        end

        rho = x * x';
    otherwise
        util.error(['Invalid equation of motion : ' eom]);
end
