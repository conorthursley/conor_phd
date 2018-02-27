% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef fft
    
    properties (Access = public)
        
        % General properties
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass that enters the kinetic energy
        periodic    % use periodic boundary conditions or not
        nokin       % enable/disable the kinetic energy operator
        
        % Grid in position space
        n_pts       % number of grid points (equal in position and momentum grid)
        x_min       % lower bound of the position grid
        x_max       % upper bound of the position grid
        
    end
    
    properties (Access = private)
        
        % Grid in position space
        x_dlt       % (constant!) grid spacing
        
        % Grid in momentum space
        p_min       % lower bound of the momentum grid
        p_max       % upper bound of the momentum grid
        
        % Internal representation of the kinetic energy operator/propagator
        kin         % DVR grid representation of the kinetic energy Hamiltonian   
        kin_shift   % an fftshifted version of kin
        kin_expo    % same, but exponentiated, for split operator
        kin_factor  % normalize and removes a spurious complex phase
        
    end
    
    methods (Access = public) % in separate files within *this* directory
        
        dvr = fbr2dvr(obj, fbr)
        fbr = dvr2fbr(obj, dvr)
        [ obj, weights, x_grid, p_grid ] = init_grid(obj)
        obj = init_kin(obj, fraction, output)
        obj = kinetic(obj, new)
        kinetic_exp(obj)
        dvrkin = kinetic2dvr(obj)
        dvr = matrix2dvr(obj, fbr)
        retval = momentum(obj, psi)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, s)
        [ trafo, weights, xgrid ] = transform(obj, n_pts, original)
        
        % Constructor method: Initialize default property values
        function obj = fft 
            obj.dof = 1; 
            obj.mass = 1; 
            obj.periodic = true; 
            obj.nokin = false; 
        end
        
        function disp (obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            util.disp ( ' ' )
            util.disp ( '**************************************************' )
            util.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
            util.disp ( 'Discretization scheme: Equally spaced grids (FFT)' )
            util.disp ( '**************************************************' )
            util.disp ( [ 'Number of grid  points : ' num2str(obj.n_pts)    ] )
            util.disp ( ' ' )
            util.disp ( 'Position space' )
            util.disp ( [ 'Minimum of grid  : ' num2str(obj.x_min)    ] )
            util.disp ( [ 'Maximum of grid  : ' num2str(obj.x_max)    ] )
            util.disp ( [ 'Grid spacing     : ' num2str(obj.x_dlt)  ] )
            util.disp ( ' ' )
            util.disp ( 'Momentum (wavenumber) space (Fourier transform)' )
            util.disp ( [ 'Maximum of grid  : ' num2str(obj.p_max)    ] )
            util.disp ( [ 'Grid spacing     : ' num2str(2*obj.p_max/obj.n_pts)  ] )
            util.disp ( ' ' )
            util.disp ( 'Momentum (wavenumber) space (Wigner transform)' )
            util.disp ( [ 'Maximum of grid  : ' num2str(obj.p_max/2)  ] )
            util.disp ( [ 'Grid spacing     : ' num2str(obj.p_max/obj.n_pts)] )
        end

    end

end
