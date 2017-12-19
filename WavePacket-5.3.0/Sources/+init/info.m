%--------------------------------------------------------------------------
%
% General information: Names of program, user, host, etc
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt
%               2007-2014 Ulf Lorenz
%
% see the README file for license details.

function info (filename)
global info

% Path and file name
[pathstr, filestr, ~] = fileparts (filename);

% Program name and version number
info.package = 'WavePacket';
info.program = filestr;

% SVN revision numbers: WavePacket installation
rev = util.revision (pathstr);
% if rev~=0
%     info.version1 = ['R' num2str(rev)];
% else
    info.version1 = 'V5.3.0';
% end

% SVN revision numbers: Current working directory
rev = util.revision ('./');
if rev~=0
    info.version2 = ['R' num2str(rev)];
else
    info.version2 = '(unversioned)';
end

% Get MATLAB release
info.release = version ('-release');

% Get user name 
result = license('inuse', 'matlab');
info.user_name = result.user;
% info.user_name = result.feature;

% Get host name
if isunix 
    [status,result] = unix('hostname');
    if ~status
        info.host_name = strcat(result,' (Unix)');
    else
        info.host_name = 'Running on Unix';
    end
elseif ispc
    [status,result] = dos('hostname');
    if ~status
        info.host_name = strcat(result,' (Windows)');
    else
        info.host_name = 'Running on Windows';
    end
end

    
% Get path and file name 
info.path_name = pwd;
addpath(info.path_name);

% Open log file (or create a new one) for writing in text mode
% Discard existing contents, if any.
filename = fullfile(info.path_name, strcat(info.program, '.log'));
fclose('all'); % needed when qm_xyz is called inside a loop
info.stdout = fopen(filename, 'wt');

if info.stdout == -1
    error(strcat('Could not open log file ',filename));
end

% Output program name/version etc
util.disp ( ' ' )
util.disp ( '************************************************************')
util.disp ( ' ' )
util.disp (['Program package : ', info.package] )
util.disp (['Program name    : ', info.program] )
util.disp (['Version number  : ', num2str(info.version1),' (MATLAB)'])
util.disp (['Version number  : ', num2str(info.version2),' (work dir.)'])
util.disp ( ' ' )
util.disp (['Path name       : ', info.path_name] )
util.disp (['User name       : ', info.user_name] )
util.disp (['Host name       : ', info.host_name] )
util.disp (['MATLAB release  : ', info.release  ] )

% Initialize stopwatch timer; Output clock/date/time
info.start_time = cputime;
util.clock;

% Output GPL and copyright info
util.disp ( '************************************************************')
util.disp ( 'This program is subject to the GNU General Public License v2')
util.disp ( 'see http://www.gnu.org or the root path for the license text')
util.disp ( '************************************************************')
util.disp ( ' ' )
util.disp ( ' (C) 2004-2017 Burkhard Schmidt, Ulf Lorenz, FU Berlin' )
util.disp ( '     http://sourceforge.net/projects/matlab.wavepacket.p' )
util.disp ( ' ' )
