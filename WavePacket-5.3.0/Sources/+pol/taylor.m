% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.

function taylor 
global hamilt space

util.disp (' ')
util.disp ('*********************************************************')
util.disp ('Polarizability: Taylor series (diagonal in N dimensions) ')
util.disp ('Optionally for x and/or y direction                      ')
util.disp ('                                                         ')
util.disp ('              inf   N   v_mn,jk              j           ')
util.disp (' alpha  (r) = Sum  Sum  ------- ( r  - s    )   + c      ')
util.disp ('      mn      j=1  k=1     j!      k    mn,k        mn   ')
util.disp ('                                                         ')
util.disp ('*********************************************************')

%% x-polarization
if isfield(hamilt.pol,'x')
    util.disp ('Polarizabilities along x-direction')
    
    if ~isfield(hamilt.pol.x,'s')
        hamilt.pol.x.s = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.pol.x,'c')
        hamilt.pol.x.c = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.pol.x,'v')
        hamilt.pol.x.v = cell (hamilt.coupling.n_eqs);
    end
    
    % so far: only diagonal elements (m=n)
    for m = 1:hamilt.coupling.n_eqs
        util.disp (['Taylor series for x-polarizability (' int2str(m) ',' int2str(m) ')'])
        
        % Reference position: default is zero
        if isempty(hamilt.pol.x.s{m,m})
            hamilt.pol.x.s{m,m} = zeros(1,space.size.n_dim); % row vector
        end
        
        % Constant offset; default is zero
        if isempty(hamilt.pol.x.c{m,m})
            hamilt.pol.x.c{m,m} = 0; % scalar
        end
        
        % Evaluate Taylor series
        if isempty (hamilt.pol.x.v{m,m})
            hamilt.p_x.grid_ND{m,m} = ...
                hamilt.pol.x.c{m,m} * ...
                ones(size(space.dvr.grid_ND{1}));
        else
            hamilt.p_x.grid_ND{m,m} = util.taylor (...
                space.dvr.grid_ND, ...
                hamilt.pol.x.s{m,m}, ...
                hamilt.pol.x.c{m,m}, ...
                hamilt.pol.x.v{m,m} );
        end
        
        util.disp ('   ')
    end
    
    else
    
    % No polarizabilities along x
    for m = 1:hamilt.coupling.n_eqs
        hamilt.p_x.grid_ND{m,m} = [];
    end

end

%% y-polarization
if isfield(hamilt.pol,'y')
    util.disp ('Polarizabilities along y-direction')
    
    if ~isfield(hamilt.pol.y,'s')
        hamilt.pol.y.s = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.pol.y,'c')
        hamilt.pol.y.c = cell (hamilt.coupling.n_eqs);
    end
    if ~isfield(hamilt.pol.y,'v')
        hamilt.pol.y.v = cell (hamilt.coupling.n_eqs);
    end
    
    % so far: only diagonal elements (m=n)
    for m = 1:hamilt.coupling.n_eqs
        util.disp (['Taylor series for y-polarizability (' int2str(m) ',' int2str(m) ')'])
        
        % Reference position: default is zero
        if isempty(hamilt.pol.y.s{m,m})
            hamilt.pol.y.s{m,m} = zeros(1,space.size.n_dim); % row vector
        end
        
        % Constant offset; default is zero
        if isempty(hamilt.pol.y.c{m,m})
            hamilt.pol.y.c{m,m} = 0; % scalar
        end
        
        % Evaluate Taylor series
        if isempty (hamilt.pol.y.v{m,m})
            hamilt.p_y.grid_ND{m,m} = ...
                hamilt.pol.y.c{m,m} * ...
                ones(size(space.dvr.grid_ND{1}));
        else
            hamilt.p_y.grid_ND{m,m} = util.taylor (...
                space.dvr.grid_ND, ...
                hamilt.pol.y.s{m,m}, ...
                hamilt.pol.y.c{m,m}, ...
                hamilt.pol.y.v{m,m} );
        end
        
        util.disp ('   ')
    end
    
else
    
    % No polarizabilities along y
    for m = 1:hamilt.coupling.n_eqs
        hamilt.p_y.grid_ND{m,m} = [];
    end

end
