% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function dvr = matrix2dvr(obj, fbr)

global space

% Very simple: Just do the correct FFT for the given degree of freedom.
% Some care has to be taken with normalisation factors, though.

if length(size(fbr)) ~= 2*space.size.n_dim
    util.error('Input is not a proper matrix.')
end


%% First step: transform the left index as in the fbr2dvr transformation.

% first, remove some factors (as in the fbr2dvr method). Unfortunately, we need it in a
% complicated matrix form, which we have to assemble and reshape and all that ourselves.
expGrid_2ND = repmat(exp(-1i * obj.x_min * space.fbr.grid_ND{obj.dof}(:)), 1, space.size.n_tot);
expGrid_2ND = reshape( expGrid_2ND, size(fbr) );
factorLeft = sqrt(obj.x_max-obj.x_min) / obj.n_pts * expGrid_2ND;

dvr = fbr ./ factorLeft;
dvr = ifft(ifftshift(dvr, obj.dof), [], obj.dof);

%% The right index is transformed as in the dvr2fbr transformation

% again, this ugly assembly of the additional phase factors.
expGrid_2ND = repmat(exp(-1i * obj.x_min * space.fbr.grid_ND{obj.dof}(:)'), space.size.n_tot, 1);
expGrid_2ND = reshape( expGrid_2ND, size(fbr) );
factorRight = sqrt(obj.x_max-obj.x_min) / obj.n_pts * expGrid_2ND;

dvr = fftshift ( fft(dvr, [], obj.dof + space.size.n_dim), obj.dof+space.size.n_dim );
dvr = dvr .* factorRight;
