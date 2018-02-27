%------------------------------------------------------------------------------
%
% This function writes a text string both to the screen and to a logfile(3)
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.
    
function disp (string)

global info

% Write text string to log file
if ~isempty(info) && isfield(info, 'stdout') && ~isempty(info.stdout)
	fprintf (info.stdout,'%s\n',string);
end

% display text string on screen
disp (string);
