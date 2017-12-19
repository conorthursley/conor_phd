% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function dvr = fbr2dvr(obj, fbr)

% Again, we have to shape the input to a pseudo-2D form before we can 
% apply the transformation matrix due to matlab limitations.
% We reshape the matrix back to the original form in the end.

[dvr, permutation, shapedims] = shape(obj, fbr);

dvr = obj.trafo2dvr * dvr;

dvr = shape_back(obj, dvr, permutation, shapedims);
