%------------------------------------------------------------------------------
%
% This function returns the (true!) number of array dimensions
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.
    
function n = ndims (array)

% Use builtin function
n = ndims(array);

% Vectors are stored as matrices
if n == 2

    % Detect column vector
    if numel(array)==size(array,1)
        n = 1;
    end

    % Detect row vector
    if numel(array)==size(array,2)
        n = 1;
    end

end
