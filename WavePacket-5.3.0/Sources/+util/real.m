%------------------------------------------------------------------------------
%
% Cast a complex value into real form unless the imaginary part becomes too
% large, then we issue a warning.
% Required sometimes, because a value that should be real by theory 
% but has a tiny (~1e-16) imaginary part due to numerical inaccuracy.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function realval = real(complexval)

realval = real(complexval);
imagval = imag(complexval);
if abs(complexval) > eps('single') && abs(imagval) > 10^-6 * abs(realval)
    util.disp(['Warning: Got a complex value where a real one was expected. Imaginary part = ' num2str(imagval)]);
end
