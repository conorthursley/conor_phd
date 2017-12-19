% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function obj = init_kin(obj, fraction, output)

global space time hamilt


if nargin < 3
    output = true;
end


obj.kin = space.fbr.grid_ND{obj.dof}.^2 / (2*obj.mass);

% Truncate the kinetic energy at some upper bound.
obj.kin(obj.kin > hamilt.truncate.delta) = hamilt.truncate.delta;

% For the propagation or application of the kinetic energy, it is much
% faster if we don't have to shift around the grids too often (internally,
% these are matrix copies!), so we fftshift the kinetic energy
% grids once here instead of applying it again and again to the wavefunction.
obj.kin_shift = ifftshift(obj.kin, obj.dof);

% Just exponentiate...
obj.kin_expo = exp(-1i * time.sub.delta * fraction * obj.kin_shift);


%% Don't know another place to put it...
% factor that we need to apply because FFT does not do it. It normalizes the
% wave function in momentum representation and removes a spurious complex phase
obj.kin_factor = sqrt(obj.x_max - obj.x_min) / obj.n_pts ...
        * exp(-1i * space.fbr.grid_ND{obj.dof} * obj.x_min);


%% Informational output

if obj.nokin
    return
end

if output
    util.disp (' ')
    util.disp ('*************************************************')
    util.disp ( ['Kinetic energy for degree of freedom: ' int2str(obj.dof)] )
    util.disp ('                                                 ')
    util.disp ('          1  (  d  ) 2                           ')
    util.disp (' T  =  - --- ( --- )                             ')
    util.disp ('         2 M ( d R )                             ')
    util.disp ('                                                 ')
    util.disp ('*************************************************')
    util.disp ( ['Mass :' num2str(obj.mass)] )
end
