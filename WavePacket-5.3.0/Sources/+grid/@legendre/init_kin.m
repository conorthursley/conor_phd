% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function obj = init_kin(obj, fraction, output)

global space time hamilt


if nargin < 3
    output = true;
end


% Builds the matrices for the kinetic operator and propagator in the 
% spectral basis.

%% First, let us build the matrix for the kinetic operator.

% The matrix is in principle diagonal in FBR, and has elements l(l+1)/2mR^2.
% In practice, this is a bit tricky, because we have to separate the cases
% "R fixed" and "R taken from another degree of freedom", so the code is 
% a tiny bit more complex.

if ~isempty(obj.R_0)
    % The easy case: rigid rotator. We directly construct the pseudo-2D form
    % here.
    obj.kin = space.fbr.grid_ND{obj.dof}.^2 / (2*obj.mass*obj.R_0^2);
else
    % The complicated case: The R value is taken from another DOF whose index
    % resides in obj.R_dof.
    obj.kin = space.fbr.grid_ND{obj.dof}.^2 ./ (2 * obj.mass * space.dvr.grid_ND{obj.R_dof}.^2);
end

%% Second: Due to all the neccessary pseudo-2D matrix reshaping, we have to
%% cast the kinetic energy operator in a pseudo-2D form as well. This can be
%% accomplished by reshaping it in the same way that we later do with the wave
%% function.

[obj.kin, a, b] = shape(obj, obj.kin);

%% Last, truncate it.
obj.kin(obj.kin > hamilt.truncate.delta) = hamilt.truncate.delta;

%% Now for the easier part: The Propagator. Just exponentiate the kinetic
%% matrix.
obj.kin_expo = exp(-1i * time.sub.delta * fraction * obj.kin);


%% Informational output

if obj.nokin
    return
end


if output
    util.disp (' ')
    util.disp ('*************************************************')
    util.disp ( ['Kinetic energy for degree of freedom: ' int2str(obj.dof)] )
    util.disp (' ')
    util.disp ('        1   -2   ^ 2                             ')
    util.disp (' T  =  --- R     L                               ')
    util.disp ('       2 M                                       ')
    util.disp (' ')
    if isempty(obj.R_dof)
        util.disp ( ['where R is chosen as constant ' num2str(obj.R_0)] )
    else
        util.disp ( ['where R is taken from the DOF ' num2str(obj.R_dof)] )
    end
    util.disp ('*************************************************')
    util.disp ( ['Mass: ' num2str(obj.mass)] )
    util.disp (' ')
end
