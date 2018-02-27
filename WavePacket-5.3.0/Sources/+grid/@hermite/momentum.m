% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = momentum(obj, psi)

% In theory this is done on every call to grid_momentum. In practise,
% the kinetic (and therefore momentum) grid should be initialised after
% the first time step as we _need_ kinetic energy for propagating.
if isempty(obj.dvr_mom)
    obj = init_kin(obj, 1);
end

% Again, shape the grid to a nice form, multiply with a matrix, and
% shape back.

[retval, permutation, shapedims] = shape(obj, psi);
retval = obj.dvr_mom * retval;
retval = shape_back(retval, permutation, shapedims);
