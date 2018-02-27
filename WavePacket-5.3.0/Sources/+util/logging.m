% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2010 Ulf Lorenz
%
% see the README file for license details.

function logging ( step )
global expect info space time uncert hamilt psi

%% Construct header from several text strings
if strcmpi(info.program,'qm_bound') % TISE: Loop over eigenstates
    info.header1 = strcat('n = ', ...
        int2str(psi.eigen.start + step-1)); % Start counting  from zero (for bound states)
else    % qm_propa or hand-crafted scripts
    info.header1 = strcat('t = ', ...
        num2str(time.main.grid(step), '%10.4f'), ...
        ', step =', num2str(step-1));
end

if hamilt.coupling.n_eqs==1
    estr = 'E = ';
else
    if hamilt.coupling.representation=='adi'
        estr = 'E_{adi} = ';
    elseif hamilt.coupling.representation=='dia'
        estr = 'E_{dia} = ';
    end
end
    
info.header2 = strcat(estr, ...
    num2str(expect.tot.tot(step), '%12.6f'), ...
    ', N = ', ...
    num2str(expect.tot.pop(step), '%10.4f'));

util.disp ('*************************************************************')
util.disp (info.header1);
util.disp (info.header2);
util.disp ('*************************************************************')
util.disp (' ')

%% Expectation values of individual wavefunctions
for m=1:hamilt.coupling.n_eqs

    % If population exceeds threshold
    if expect.ind.pop{m}(step)>expect.min_pop

        % Population
        if hamilt.coupling.n_eqs>1
            if strcmpi(hamilt.coupling.representation,'adi')
                util.disp (['Adiabatic state ',int2str(m),': Population ',num2str(expect.ind.pop{m}(step),'%10.4f')])
            elseif strcmpi(hamilt.coupling.representation,'dia')
                util.disp (['Diabatic state ',int2str(m),' (',hamilt.coupling.labels{m},') : Population ',num2str(expect.ind.pop{m}(step),'%10.4f')])
            end
            util.disp ('-------------------------------------------------')
            util.disp (' ')
        end

        % additional multiplicative operators
        if isfield(space, 'amo') 
            for p=1:length(space.amo)
                util.disp ([space.amo{p}.label ': ',num2str(expect.ind.amo{m}(step,p),'%12.6f')])
            end
            util.disp (' ')
        end
        
        % Position/momentum
        util.disp (['Position ',num2str(expect.ind.dvr{m}(step,:),'%12.6f'),' +/- ',num2str(uncert.ind.dvr{m}(step,:),'%12.6f')])
        util.disp (['Momentum ',num2str(expect.ind.fbr{m}(step,:),'%12.6f'),' +/- ',num2str(uncert.ind.fbr{m}(step,:),'%12.6f')])
        util.disp (['Uncertainty product ',num2str(uncert.ind.dvr{m}(step,:).*uncert.ind.fbr{m}(step,:),'%12.6f')])
        util.disp (' ')

        % Potential/kinetic/sum energy
        util.disp (['Potential energy ',num2str(expect.ind.pot{m}(step),'%15.8f'),' +/- ',num2str(uncert.ind.pot{m}(step),'%12.6f')])
        util.disp (['Kinetic energy   ',num2str(expect.ind.kin{m}(step),'%15.8f'),' +/- ',num2str(uncert.ind.kin{m}(step),'%12.6f')])
        util.disp (['Sum of energies  ',num2str(expect.ind.all{m}(step),'%15.8f')])
        util.disp (' ')

    end
end

%% Expectation values of total  wavefunctions
if hamilt.coupling.n_eqs>1

    % Population
    util.disp (['Total: Population ',num2str(expect.tot.pop(step),'%10.4f')])
    util.disp ('-------------------------------------------------------------')
    util.disp (' ')

    % Projection
    if isfield(space, 'amo') 
        for p=1:length(space.amo)
            util.disp ([space.amo{p}.label ': ',num2str(expect.tot.amo(step,p),'%12.6f')])
        end
        util.disp (' ')
    end

    % Position/momentum
    util.disp (['Position ',num2str(expect.tot.dvr(step,:),'%12.6f')])
    util.disp (['Momentum ',num2str(expect.tot.fbr(step,:),'%12.6f')])
    util.disp (' ')

    % Potential/kinetic/sum energy
    util.disp (['Potential energy ',num2str(expect.tot.pot(step),'%12.6f')])
    util.disp (['Kinetic energy   ',num2str(expect.tot.kin(step),'%12.6f')])
    util.disp (['Sum of energies  ',num2str(expect.tot.all(step),'%12.6f')])
    util.disp (' ')
    util.disp (['Incl. couplings  ',num2str(expect.tot.tot(step),'%12.6f')])
    util.disp (' ')

end

