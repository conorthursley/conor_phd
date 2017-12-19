%--------------------------------------------------------------------------
% Propagate wavefunction (psi.dvr.grid_ND{m}) by time.sub.n substeps of size 
% time.sub.delta subject to a given Hamiltonian (hamilt) by repeated 
% action of finite differences
%
%                            i      ^              i    ^3
% psi(t+dt) = psi(t-dt) - 2 ---- dt H psi(t) + -------- H  psi(t) - ...
%                           hbar               3 hbar^3
%
% Only second order differencing is implemented up to now, which has a
% fairly poor convergence.
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

function differencing ( step )
global psi space time hamilt

% First step only: Initialize "old" wavefunction; no propagation yet
if step==1
    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ('Numerical propagator: Finite differences')
    util.disp ('*******************************************************')

    if ~isfield (time.propa,'order')
        time.propa.order = 2;
    end
    
    switch time.propa.order
        case (2)
            util.disp ( 'Second order differencing (SOD)' )
        otherwise
            util.error ( 'Wrong choice of differencing order' )
    end
    util.disp ( ' ' )
    psi.dvr.old_ND = psi.dvr.grid_ND;
else
    
    % All later steps: Perform propagation
    % Loop over substeps
    for k = 1:time.sub.n

        % Perform propagation for one substep
        if time.efield.n_pulse>0 
            ket.hamilt(time.efield.short.x(k), time.efield.short.y(k), 0);
        else
            ket.hamilt(0,                      0,                      0);
        end
        
        for m = 1:hamilt.coupling.n_eqs
            psi.dvr.new_ND{m} = psi.dvr.old_ND{m} - 2 * 1i * time.sub.delta * psi.dvr.new_ND{m};
        end
            
        % Get ready for next time step
        psi.dvr.old_ND = psi.dvr.grid_ND;
        psi.dvr.grid_ND = psi.dvr.new_ND;

        % Autocorrelation function (using position representation)
        time.acf.short(k) = 0;
        for m = 1:hamilt.coupling.n_eqs
            time.acf.short(k) = time.acf.short(k) + ...
                sum ( conj(psi.dvr.init.ND{m}(:)).*psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:) );
        end
                
    end

end

