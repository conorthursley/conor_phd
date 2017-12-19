% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function error (string)

% Write error message to log file
fprintf (3,'%s\n',string);

% Display error message on screen and terminate
error (string);
