% Create input signal for simulation:
input = frest.Sinestream('Frequency',logspace(1,11,40000));

% Open the Simulink model:
open TwoCellAMM

% Specify portion of model to estimate:
io(1)=linio('TwoCellAMM/Constant',1,'input');
io(2)=linio('TwoCellAMM/Transfer Fcn2',1,'openoutput');

% Specify the steady state operating point for the estimation.
% TwoCellAMM_spec = operspec('TwoCellAMM');
% op = findop('TwoCellAMM',TwoCellAMM_spec);

% Estimate frequency response of specified blocks:
sysest = frestimate('TwoCellAMM',io,input);
bode(sysest)