% Kinetic energy in Jacobi coordinates for the angle.
% See for example J. Chem. Phys 116:4403
 
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef jacobi
       
    properties (Access = public)
        
        r_0   
        dof_R 
        dof_r 
        dof_c        
        mass_R
        mass_r
        
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
