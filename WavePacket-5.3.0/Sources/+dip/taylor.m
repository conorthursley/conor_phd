% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 Burkhard Schmidt
%
% see the README file for license details.

function taylor 
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Dipole moments: Taylor series (diag. in N dimensions)')
util.disp ('Optionally for x and/or y direction                  ')
util.disp ('                                                     ')
util.disp ('          inf   N   v_mn,jk              j           ')
util.disp (' d  (r) = Sum  Sum  ------- ( r  - s    )   + c      ')
util.disp ('  mn      j=1  k=1     j!      k    mn,k        mn   ')
util.disp ('                                                     ')
util.disp ('*****************************************************')

%% x-polarization
if isfield(hamilt.dip,'x')
    util.disp ('Dipole moments along x-direction')
    
    if ~isfield(hamilt.dip.x,'s')
        hamilt.dip.x.s = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.dip.x,'c')
        hamilt.dip.x.c = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.dip.x,'v')
        hamilt.dip.x.v = cell (hamilt.coupling.n_eqs);
    end
    
    % permanent (m=n) and transition dipole moments ( m<n )
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            util.disp (['Taylor series for x-dipole(' int2str(m) ',' int2str(n) ')'])
            
            % Reference position: default is zero
            if isempty(hamilt.dip.x.s{m,n})
                hamilt.dip.x.s{m,n} = zeros(1,space.size.n_dim); % row vector
            end
            
            % Constant dipole offset; default is zero
            if isempty(hamilt.dip.x.c{m,n})
                hamilt.dip.x.c{m,n} = 0; % scalar
            end
            
            % Evaluate Taylor series only for desired {m,n} combinations
            if isempty (hamilt.dip.x.v{m,n})
                hamilt.d_x.grid_ND{m,n} = ...
                    hamilt.dip.x.c{m,n} * ...
                    ones(size(space.dvr.grid_ND{1}));
            else
                hamilt.d_x.grid_ND{m,n} = util.taylor (...
                    space.dvr.grid_ND, ...
                    hamilt.dip.x.s{m,n}, ...
                    hamilt.dip.x.c{m,n}, ...
                    hamilt.dip.x.v{m,n} );
            end
            
            util.disp ('   ')
        end
    end
    
else
    
    % No dipole moments along x
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            hamilt.d_x.grid_ND{m,n} = [];
        end
    end
   
end

%% y-polarization
if isfield(hamilt.dip,'y')
    util.disp ('Dipole moments along y-direction')
    
    if ~isfield(hamilt.dip.y,'s')
        hamilt.dip.y.s = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.dip.y,'c')
        hamilt.dip.y.c = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.dip.y,'v')
        hamilt.dip.y.v = cell (hamilt.coupling.n_eqs);
    end
    
    % permanent (m=n) and transition dipole moments ( m<n )
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            util.disp (['Taylor series for y-dipole(' int2str(m) ',' int2str(n) ')'])
            
            % Reference position: default is zero
            if isempty(hamilt.dip.y.s{m,n})
                hamilt.dip.y.s{m,n} = zeros(1,space.size.n_dim); % row vector
            end
            
            % Constant dipole offset; default is zero
            if isempty(hamilt.dip.y.c{m,n})
                hamilt.dip.y.c{m,n} = 0; % scalar
            end
            
            % Evaluate Taylor series only for desired {m,n} combinations
            if isempty (hamilt.dip.y.v{m,n})
                hamilt.d_y.grid_ND{m,n} = ...
                    hamilt.dip.y.c{m,n} * ...
                    ones(size(space.dvr.grid_ND{1}));
            else
                hamilt.d_y.grid_ND{m,n} = util.taylor (...
                    space.dvr.grid_ND, ...
                    hamilt.dip.y.s{m,n}, ...
                    hamilt.dip.y.c{m,n}, ...
                    hamilt.dip.y.v{m,n} );
            end
            
            util.disp ('   ')
        end
    end
    
else
    
    % No dipole moments along y
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            hamilt.d_y.grid_ND{m,n} = [];
        end
    end    
    
end
