%--------------------------------------------------------------------------
%
% This routine reconstructs a wavefunction from a plane wave expansion
% using FFTW as included in Matlab) to transform from FBR to DVR
% representation. See grid_dvr2fbr for some more extensive documentation
% of what we do.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function dvr = fbr2dvr(obj, fbr)

% Note that we have introduced some prefactors in the expansion, which we have
% to get rid of again.

dvr = fbr ./ obj.kin_factor;
dvr = ifft  ( ifftshift ( dvr, obj.dof ), [], obj.dof );
