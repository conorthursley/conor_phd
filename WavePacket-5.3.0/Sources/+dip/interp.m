%--------------------------------------------------------------------------
% 
% Read tablulated values of dipole moment functions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation. Using 
% separate files for x and y components!
% 
% Entries of dipole moment matrices:
% Permanent dipole moments: 'd_x_m.dat', 'd_y_m.dat' with m=1,2,...
% Transition dipole moments: 'd_x_m_n.dat', 'd_y_m_n.dat' with n>m
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
util.disp ('*******************************************************')
util.disp ('Interpolation of tabulated values of dipole moments')
util.disp ('*******************************************************')

if ~isfield (hamilt.dip,'pos_conv')
    hamilt.dip.pos_conv = 1;
end
util.disp (['Conversion factor for coordinates    : ' num2str(hamilt.dip.pos_conv)])

if ~isfield (hamilt.dip,'dip_conv')
    hamilt.dip.dip_conv = 1;
end
util.disp (['Conversion factor for dipole moments : ' num2str(hamilt.dip.dip_conv)])

if ~isfield (hamilt.dip,'method')
    hamilt.dip.method = 'spline';
end
util.disp (['Interpolation method                 : '         hamilt.dip.method])

if space.size.n_dim>1
    if length(hamilt.dip.n_pts)~=space.size.n_dim
        error ('Inconsistent number of dimensions of tabulated data')
    end
    util.disp (['Number of data points expected    : ' int2str(hamilt.dip.n_pts)])
else
    hamilt.dip.n_pts = 1; % Dummy
end
util.disp ( ' ' )

% Triple loop over entries of dipole moment matrix (x,y components)
for m=1:hamilt.coupling.n_eqs
    for n=m:hamilt.coupling.n_eqs
        for p=['x' 'y']

            % Name of input data file
            if m==n
                filename = strcat('d_', p, '_', int2str(m), '.dat');
                util.disp ([p, '-Dipole : ' hamilt.coupling.labels{m}])
            else
                filename = strcat('d_', p, '_', int2str(m), '_', int2str(n), '.dat');
                util.disp ([p, '-Dipole: ' hamilt.coupling.labels{m},' <--> ',hamilt.coupling.labels{n}])
            end

            % Check availability of input data file
			fid = fopen(filename);
            if fid == -1

                util.disp (['Input data file not available : ' filename])
                util.disp ('Leaving entry of dipole matrix empty')
                util.disp (' ')
                if p=='x'
                    hamilt.d_x.grid_ND{m,n} = [];
                elseif p=='y'
                    hamilt.d_y.grid_ND{m,n} = [];
                end
            else
				fclose(fid);

                if p=='x'
                    hamilt.d_x.grid_ND{m,n} = util.interp (...
                        space.dvr.grid_ND, ...
                        filename, ...
                        hamilt.dip.n_pts, ...
                        hamilt.dip.method, ...
                        hamilt.dip.pos_conv, ...
                        hamilt.dip.dip_conv);
                elseif p=='y'
                    hamilt.d_y.grid_ND{m,n} = util.interp (...
                        space.dvr.grid_ND, ...
                        filename, ...
                        hamilt.dip.n_pts, ...
                        hamilt.dip.method, ...
                        hamilt.dip.pos_conv, ...
                        hamilt.dip.dip_conv);
                end

            end

            util.disp(' ')

        end
    end
end
