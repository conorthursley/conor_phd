%--------------------------------------------------------------------------
% Visualize spectrum, wavefunctions, and potential 
% from previous qm_bound simulations
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 by Burkhard Schmidt
%
% see the README file for license details.

function qm_plot(savedir, savefile)

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Provide default values for missing input arguments
if nargin<2
    savefile = 'WavePacketSave';
end

if nargin<1
    savedir=pwd;
end

util.disp('-------------------------------------------------');
util.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_movie')
util.disp('Plot wavefunctions from saved calculation in file        ');
util.disp(['     ' savefile]);
util.disp('residing in directory                            ');
util.disp(['     ' savedir]);
util.disp('-------------------------------------------------');
util.disp(' ');


figure(88);
clf



    S = load(filename);
    T = load(strcat(filename,'_0'));
    
    % Plot potential energy function
    r = 
    pot = S.hamilt.pot.grid_ND{1,1};
    plot (r, pot,'k')
    hold on
    
    % Eigen-energies and -functions
    for n=S.psi.eigen.start:S.psi.eigen.stop
        ene = S.hamilt.eigen.eig_vals(n+1);
        psi = T.save_wf{n+1}{1};
        
        
        % Plot lines/curves
            plot(theta, ene*ones(size(theta)),'b--')
            plot(theta, ene+psi*5,'b')
        
    end
    hold off
    
	
	% Output clock/date/time
util.clock;

end
 


