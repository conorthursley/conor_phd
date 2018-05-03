%% Script to run the other scripts in this folder over a weekend
% each script takes about 24 hrs to run and I'm doing 3 of them so the need
% for a weekend.

% each script should have the results and figure saved. 
% the figure is saved because all of the useful data is in that figure,
% which can be extracted for later use.
close all
clear all

%% Run 2U system code with error handling capabilities
try 
    NL_MIM_2U_SweptSine    % no need to add the .m at the end of it
    figureFile2='U:\_PhD\Matlab\Git\Models\MDOF Metamaterial Models\SweptSine\WeekendSimulationCode\2UweaklyNL.fig';
    sendolmail('conor.macdonald@adelaide.edu.au','2U simulation complete',...
    '2U simulation with figure',{figureFile2});
catch 
    warning('Problem using the 2U simulation. Sending email.')
    exception = MException.last;
    msgtext = getReport(exception);
    sendolmail('conor.macdonald@adelaide.edu.au','2U simulation failed',...
    ['2U failed due to ' msgtext newline ' -----------:End message:--------']);
end
%% Run 5U system code

try 
    NL_MIM_5U_SweptSine    % no need to add the .m at the end of it
    figureFile5='U:\_PhD\Matlab\Git\Models\MDOF Metamaterial Models\SweptSine\WeekendSimulationCode\5UweaklyNL.fig';
    sendolmail('conor.macdonald@adelaide.edu.au','5U simulation complete',...
    '5U simulation with figure',{figureFile5});
catch 
    warning('Problem using the 5U simulation. Sending email.')
    exception = MException.last;
    msgtext = getReport(exception);
    sendolmail('conor.macdonald@adelaide.edu.au','5U simulation failed',...
    ['5U failed due to ' msgtext newline ' -----------:End message:--------']);
end
%% Run 10U system code

try 
    NL_MIM_10U_SweptSine    % no need to add the .m at the end of it
    figureFile10='U:\_PhD\Matlab\Git\Models\MDOF Metamaterial Models\SweptSine\WeekendSimulationCode\10UweaklyNL.fig';
    sendolmail('conor.macdonald@adelaide.edu.au','10U simulation complete',...
    '10U simulation with figure',{figureFile10});
catch 
    warning('Problem using the 10U simulation. Sending email.')
    exception = MException.last;
    msgtext = getReport(exception);
    sendolmail('conor.macdonald@adelaide.edu.au','10U simulation failed',...
    ['10U failed due to ' msgtext newline ' -----------:End message:--------']);
end
