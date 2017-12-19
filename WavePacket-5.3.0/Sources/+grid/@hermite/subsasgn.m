% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = subsasgn(obj, index, val)

global space

if index.type == '.'
    switch index.subs
        case 'label'
            obj.label = val;
        case 'dof'
            if val > space.size.n_dim
                util.error ('too large value for index "dof" assigned');
            end
			obj.dof = val;
            % If necessary, convert obj.dof into a label for this d.o.f.
            if isempty(obj.label)
                if space.size.n_dim > 1
                    obj.label = int2str(obj.dof);
                else
                    obj.label = '';
                end
            end
		case 'mass'
			if val <= 0
				util.error('non-positive mass assigned');
			end
			obj.mass = val;
        case 'omega'
            if val <= 0
                util.error('non-positive frequency assigned');
            end
            obj.omega = val;
            obj.v_2 = [];
        case 'v_2'
            if val <= 0
                util.error('non-positive force constant assigned');
            end
            obj.v_2 = val;
            obj.omega = [];
        case 'r_e'
            obj.r_e = val;
        case 'n_pts'
            obj.n_pts = val;
        case 'nokin'
            if ~islogical(val) || length(val) ~= 1
                util.err('Value has to be set to a single boolean.')
            end
            obj.nokin = val;
		otherwise
			util.error( ['Tried to set invalid field ' index.subs ' of Hermite DVR'] );
	end
end
