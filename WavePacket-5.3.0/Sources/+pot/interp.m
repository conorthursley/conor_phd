%--------------------------------------------------------------------------
% 
% Read tabulated values of potential energy functions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation 
% 
% Entries of diabatic potential energy matrices:
% Diagonal elements: 'pot_m.dat' with m=1,2,...
% Off-diagonal elements: 'pot_m_n.dat' with n>m
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function interp
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Interpolation of tabulated values of potential energy')
util.disp ('*****************************************************')

if ~isfield (hamilt.pot,'pos_conv')
    hamilt.pot.pos_conv = 1;
end
util.disp (['Conversion factor for coordinates : ' num2str(hamilt.pot.pos_conv)])

if ~isfield (hamilt.pot,'pot_conv')
    hamilt.pot.pot_conv = 1;
end
util.disp (['Conversion factor for energies    : ' num2str(hamilt.pot.pot_conv)])

if ~isfield (hamilt.pot,'method')
    hamilt.pot.method = 'spline';
end
util.disp (['Interpolation method              : '         hamilt.pot.method])

if space.size.n_dim>1
    if length(hamilt.pot.n_pts)~=space.size.n_dim
        error ('Inconsistent number of dimensions of tabulated data')
    end
    util.disp (['Number of data points in each dimension  : ' int2str(hamilt.pot.n_pts)])
else
    hamilt.pot.n_pts = 1; % Dummy
end
util.disp ( ' ' )

% Double loop over entries of (diabatic) potential energy matrix
for m=1:hamilt.coupling.n_eqs
    for n=m:hamilt.coupling.n_eqs
        
        % Name of input data file
        if m==n
            filename = strcat('pot', '_', int2str(m), '.dat');
            if hamilt.coupling.n_eqs>1
                util.disp (['Potential: ' hamilt.coupling.labels{m}])
            end
        else
            filename = strcat('pot', '_', int2str(m), '_', int2str(n), '.dat');
            util.disp (['Potential: ' hamilt.coupling.labels{m},' <--> ',hamilt.coupling.labels{n}])
        end

        % Check availability of input data file
		fid = fopen(filename);
        if fid == -1
            if m==n
                util.error (['Input data file missing : ' filename])
            else
                util.disp (['Input data file not available : ' filename])
                util.disp ('Leaving entry of potential matrix empty')
                util.disp (' ')
                hamilt.pot.grid_ND{m,n} = [];
            end
        else
			fclose(fid);
            hamilt.pot.grid_ND{m,n} = util.interp (...
                space.dvr.grid_ND, ...
                filename, ...
                hamilt.pot.n_pts, ...
                hamilt.pot.method, ...
                hamilt.pot.pos_conv, ...
                hamilt.pot.pot_conv);
        end 
           
        util.disp(' ')
        
    end % for n
end % for m
