% The Hermite DVR is suited for an expansion in eigenfunctions of the harmonic
% oscillator (Hermite polynomials times an exponential). In contrast to other DVR
% methods, the grid points are not natively bounded to some interval, but can lie
% anywhere (of course, they are always located somewhere around the minimum of the
% harmonic oscillator).
% 
% Furthermore, the Gauss-Hermite quadrature is conveniently defined in scaled
% coordinates (corresponding to m * omega = 1, r_e = 0).  To use them for real
% problems, you have to supply the shift and the properties of the harmonic
% potential you want to use.

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef hermite
    
    properties (Access = public)
        
        % General properties
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass that enters the kinetic energy
        nokin       % enable/disable the kinetic energy operator
        
        % Harmonic oscillator properties
        omega       % angular frequency of the potential: k = m * omega^2.
        v_2         % force constant
        
        % Grid options
        n_pts       % number of points (== number of basis functions)
        r_e         % equilibrium position of the harmonic oscillator.
        
    end
    
    properties (Access = private)
        
        % Internal representation of the kinetic energy operator/propagator
        kin         % DVR grid representation of the kinetic energy Hamiltonian
        kin_expo    % Same but exponentiated for split operator
        kin_max 
        dvr_mom     % DVR(!) grid representation of the momentum operator
        
        % Arrays for transforming the wave function
        trafo2fbr   % transformation matrix : DVR=>FBR
        trafo2dvr   % transformation matrix : FBR=>DVR
        
    end
    
    methods (Access = public) % in separate files within *this* directory
        
        fbr = dvr2fbr(obj, dvr)
        dvr = fbr2dvr(obj, fbr)
        [ obj, weights, x_grid, p_grid ] = init_grid(obj)
        obj = init_kin(obj, fraction, output)
        obj = kinetic(obj, usenew)
        kinetic_exp(obj)
        dvrkin = kinetic2dvr(obj)
        dvr = matrix2dvr(obj, fbr)
        retval = momentum(obj, psi)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, index)
        [ trafo, weights, xgrid ] = transform(obj, n_pts, original)
        
        % Constructor method: Initialize default property values
        function obj = hermite
            obj.dof   = 1;
            obj.mass  = 1;
            obj.omega = 1;
            obj.r_e   = 0;
            obj.nokin = false;
            obj.kin_max = 0;
        end
        
        function disp (obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            util.disp ( ' ' )
            util.disp ( '**************************************************' )
            util.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
            util.disp ( 'Discretization scheme: Gauss-Hermite DVR' )
            util.disp ( '**************************************************' )
            util.disp ( ' ' )
            util.disp ( 'Parameters of the grid' )
            util.disp ( ['Number of grid points : ' int2str(obj.n_pts)])
            util.disp ( ['scaling factor        : ' num2str(sqrt(obj.mass*obj.omega))])
            util.disp ( ['harmonic frequency    : ' num2str(obj.omega)])
            util.disp ( ['position of minimum   : ' num2str(obj.r_e)  ])
        end

        
    end
    
    methods (Access = private) % in separate files within *private* directory
        
        [weights, xi] = quadrature (num_points)
        [shapedmat, permutation, shapedims] = shape(obj, input)
        reshaped_mat = shape_back(input, permutation, shapedims)
                
    end
    
end
