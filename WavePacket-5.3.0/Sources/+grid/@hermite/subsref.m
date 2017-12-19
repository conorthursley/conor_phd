% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, index)

global space

if index.type == '.'
	switch index.subs
        % first public properties that are required from all grid classes
		case 'label'
			retval = obj.label;
		case 'dof'
			retval = obj.dof;
        case 'n_pts'
            retval = obj.n_pts;
        case 'dvr_min'
            retval = min(space.dvr.grid_1D{obj.dof}(:));
        case 'dvr_max'
            retval = max(space.dvr.grid_1D{obj.dof}(:));
        case 'fbr_min'
            retval = 0;
        case 'fbr_max'
            retval = obj.n_pts - 1;
        case 'kin_max'
            if obj.nokin
                retval = 0;
            elseif obj.kin_max > 0
                retval = obj.kin_max;
            else
                % energy is omega * (n+0.5), which is not quite correct, since
                % we only use the_kinetic_ energy operator internally, but, well...
                % in lack of better numbers
                retval = obj.omega * (obj.n_pts - 0.5);
            end
        case 'nokin'
            retval = obj.nokin;


        % Next the Hermite-specific properties
        case 'mass'
            retval = obj.mass;
        case 'omega'
            retval = obj.omega;
        case 'v_2'
            retval = obj.v_2;
        case 'r_e'
            retval = obj.r_e;

        % Finally the private properties
        case 'dvr_mom'
            retval = obj.dvr_mom;
		case 'kin'
			retval = obj.kin;
		case 'kin_expo'
			retval = obj.kin_expo;
        case 'trafo2fbr'
            retval = obj.trafo2fbr;
        case 'trafo2dvr'
            retval = obj.trafo2dvr;
		otherwise
			util.error( ['unknown field ' index.subs ' of grid.hermite requested']);
	end
else
	util.error('Tried unsupported access form for class grid.hermite');
end
