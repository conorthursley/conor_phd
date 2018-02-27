% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function obj = subsasgn(obj, index, val)
    
global space

% We allow only setting of the parameters that are needed,
% otherwise always return an error.
if index.type == '.'
    switch index.subs
        case 'dof_R'
            obj.dof_R = val;
        case 'dof_r'
            obj.r_0 = [];       % make sure only one thing is set.
            obj.dof_r = val;
        case 'dof_c'
            obj.dof_c = val;
            space.dof{obj.dof_c}.nokin = true;
        case 'r_0'
            obj.dof_r = [];
            obj.r_0 = val;
        case 'mass_R'
            obj.mass_R = val;
        case 'mass_r'
            obj.mass_r = val;
        otherwise
            util.error ( ['Tried to set invalid field ' index.subs 'in kinetic object'] );
    end
end
