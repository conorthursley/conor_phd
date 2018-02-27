% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef triatomic
    
    properties (Access = public)
        
        dof 
        mass
        theta
        
    end
    
    properties (Access = private)
        
        grid
        grid_exp
        
    end
    
    methods (Access = public)
        
        obj = apply(obj, new)
        apply_exp(obj)
        obj = init_kin(obj, fraction, output)
        dvrkin = kinetic2dvr(obj)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, s)
        
    end

end
