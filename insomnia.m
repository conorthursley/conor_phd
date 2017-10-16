% insomnia
% =========================
%
% Little utility function to prevent the computer from entering sleep mode
% on Windows systems, while doing MATLAB computations!
% Inspired by the original Insomnia tool written in C .NET :
% http://blogs.msdn.com/b/delay/archive/2009/09/30/give-your-computer-insomnia-free-tool-and-source-code-to-temporarily-prevent-a-machine-from-going-to-sleep.aspx
%
% Example:
%
%   insomnia('on','verbose');     % force the computer to remain awake
%     ... do long computation which should not be interrupted by sleep mode...
%   insomnia('off','verbose');    % allow the computer enter sleep mode if necessary
%
% by Francesco Montorsi, 15/04/2012
%
function insomnia(enter_sleep,verbose)
    persistent insomnia_prevExecState
    
    % enter_sleep may be a string 'on' or 'off' or directly a boolean:
    if ischar(enter_sleep)
        if strcmpi(enter_sleep,'on')
            enter_sleep = true;
        else
            enter_sleep = false;
        end
    end
        
    % verbose may be a string 'verbose' or 'quiet' or directly a boolean:
    if ischar(verbose)
        if strcmpi(verbose,'verbose')
            verbose = true;
        else
            verbose = false;
        end
    end
    
    if ispc
        if enter_sleep
            % see http://msdn.microsoft.com/en-us/library/windows/desktop/aa373208(v=vs.85).aspx
            execstate.ES_CONTINUOUS = hex2dec('80000000');
            execstate.ES_SYSTEM_REQUIRED = hex2dec('00000001');

            % we need to tell MATLAB where to find insomnia_header.h
            thisfile = mfilename('fullpath');
            [pathstr, a, b] = fileparts(thisfile);
            header_fname = [pathstr '\insomnia_header.h'];
            
            % we need to call a function from kernel32.dll:
            addpath([getenv('SYSTEMROOT') '\system32']);
            if ~libisloaded('kernel32')
                % load the kernel32.dll
                loadlibrary('kernel32',header_fname);
                if ~libisloaded('kernel32')
                    fprintf('Error while trying to load kernel32.dll\n');
                    return;
                end
            end
            
            % prevent sleep mode
            insomnia_prevExecState = calllib('kernel32','SetThreadExecutionState', ...
                                    execstate.ES_CONTINUOUS | execstate.ES_SYSTEM_REQUIRED);
            if insomnia_prevExecState == 0
                fprint('Error while trying to prevent system sleep!\n');
                return;
            end
            
            if nargin>1 && verbose
                fprintf('Successfully forced computer awake mode.\n');
            end
        else
            if ~libisloaded('kernel32')
                fprintf('kernel32.dll was unloaded?\n');
                return;
            end
            
            % restore previous
            ret = calllib('kernel32','SetThreadExecutionState', insomnia_prevExecState);
            if ret == 0
                fprint('Error while trying to restore thread execution state!\n');
                return;
            end
            
            if nargin>1 && verbose
                fprintf('Successfully restored computer power policy.\n');
            end
        end
    else
        % TODO?
        fprintf('Insomnia has no implementation for this platform :(\n');
    end
end