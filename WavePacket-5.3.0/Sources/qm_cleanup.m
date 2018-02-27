%------------------------------------------------------------------------------
%
% Cleans up the calculation things (mainly closes log file for now).
%
% https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_cleanup
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% see the README file for license details.

function qm_cleanup()

global info

fclose (info.stdout);
info.stdout = [];

rmpath(info.path_name);
