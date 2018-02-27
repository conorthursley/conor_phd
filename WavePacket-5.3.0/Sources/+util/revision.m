%----------------------------------------------------------------
% Get SVN revision number for a given path 
%
% inspired by Kurt Fischle's function SubWCRev.m
% see Matlab file exchange #47520 (2014/15)
%
%----------------------------------------------------------------
function number = revision (pathstr)

% If *not* run on a Windows (PC) platform, return a zero
% Code for Unix/Mac to be added later
if ~ispc
    number = 0;
    return
end

% Try executing SubWCRev executable which comes with TortoiseSVN
% If SubWCRev doesn't exist (or fails otherwise), return a zero 
command = ['SubWCRev.exe ' pathstr];
[status,cmdout] = system (command);
if status ~= 0
    number = 0;
    return
end
 
% Break output into lines and retrieve
% revision number of last commit from second line 
lines = regexp( cmdout, '\n', 'split' );
last = regexp(lines{2},'Last committed at revision (\d+)','tokens');
number = last{1}{1};