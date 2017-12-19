% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, s)

global space

% For now, allow access to all internal data
if s.type == '.'
    switch s.subs
        case 'dof'
            retval = obj.dof;
        case 'mass'
            retval = obj.mass;
        case 'theta'
            retval = obj.theta;

        % required: upper bound for kinetic energy
        case 'kin_max'
            mass_A = obj.mass(1);
            mass_B = obj.mass(2);
            mass_C = obj.mass(3);
            mass_AB = 1 / ( 1/mass_A + 1/mass_B );
            mass_BC = 1 / ( 1/mass_B + 1/mass_C );

            retval = space.dof{obj.dof(1)}.fbr_max^2 / (2*mass_AB) ...
                    +space.dof{obj.dof(2)}.fbr_max^2 / (2*mass_BC) ...
                    +cos(obj.theta) / (2*mass_B) * space.dof{obj.dof(1)}.fbr_max ...
                        * space.dof{obj.dof(2)}.fbr_max;

        % Also return the grids, though this should never be used
        case 'grid'
            retval = obj.grid;
        case 'grid_exp'
            retval = obj.grid_exp;
        otherwise
            util.error ( ['Tried to get invalid field ' s.subs 'in kinetic object'] );
    end
else
    util.error( 'Tried unsupported access form of class kin.triatomic' );
end
