% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, s)

global space

% For now, allow access to all internal data
if s.type == '.'
    switch s.subs
        case 'dof_R'
            retval = obj.dof_R;
        case 'dof_r'
            retval = obj.dof_r;
        case 'dof_c'
            retval = obj.dof_c;
        case 'r_0'
            retval = obj.r_0;
        case 'mass_R'
            retval = obj.mass_R;
        case 'mass_r'
            retval = obj.mass_r;

        % Returns the worst case maximum kinetic energy
        case 'kin_max'
            retval = 1 / ( 2*obj.mass_R*min(space.dvr.grid_1D{obj.dof_R}(:))^2 );

            if isempty(obj.r_0)
                retval = retval + 1 / ( 2*obj.mass_r*min(space.dvr.grid_1D{obj.dof_r}(:))^2 );
            else
                retval = retval + 1 / ( 2*obj.mass_r*obj.r_0^2 );
            end

            lmax = space.dof{obj.dof_c}.l_max;
            retval = retval * lmax * (lmax + 1);    

        % Also return the grids, though this should never be used
        % except for debugging and similar things.
        case 'grid'
            retval = obj.grid;
        case 'grid_exp'
            retval = obj.grid_exp;
        otherwise
            util.error ( ['Tried to get invalid field ' s.subs 'in kin.jacobi object'] );
    end
else
    util.error( 'Tried unsupported access form of class kin.jacobi' );
end
