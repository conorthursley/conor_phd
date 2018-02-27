% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = subsasgn(obj, index, val)
    
% We allow only setting of the parameters that are needed,
% otherwise always return an error.
if index.type == '.'
    switch index.subs
        case 'dof'
            if length(val) ~= 2
                util.error('Operator needs to operate on two dimensions');
            end
            obj.dof = val;
        case 'mass'
            if length(val) ~= 3
                util.error('Operator needs three masses for setup');
            end
            obj.mass = val;
        case 'theta'
            obj.theta = val;
        otherwise
            util.error ( ['Tried to set invalid field ' index.subs 'in kinetic object'] );
    end
end
