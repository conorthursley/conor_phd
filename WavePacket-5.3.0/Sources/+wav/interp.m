%--------------------------------------------------------------------------
% 
% Read tablulated values of wavefunctions from formatted data 
% file(s) and perform (one- or multi-dimensional) interpolation 
% Reading real and imaginary parts from two separate files
% (both of which may be omitted in which case zeros are used)
% 
% Entries for the different states:
% wav_m.dat, wav_m_im.dat for the m-th state
%
% Missing wavefunctions are ignored.
%
% Large parts of the code are taken from pot_interp.m
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2009 Ulf Lorenz
%
% see the README file for license details.

function interp

global hamilt space psi

util.disp (' ')
util.disp ('**********************************************************')
util.disp ('Interpolation of initial wave function from tabulated data')
util.disp ('**********************************************************')

if ~isfield (psi,'corr')
    psi.corr = [];
end

if ~isfield (psi.corr,'conv')
    psi.corr.conv = 1;
end
util.disp (['Conversion factor for coordinates : ' num2str(psi.corr.conv)])

if ~isfield (psi.corr,'method')
    psi.corr.method = 'spline';
end
util.disp (['Interpolation method              : ' psi.corr.method])

if space.size.n_dim>1
    if length(psi.corr.n_pts)~=space.size.n_dim
        error ('Inconsistent number of dimensions of tabulated data')
    end
    util.disp (['Number of data points in each dimension  : ' int2str(psi.corr.n_pts)])
else
    psi.corr.n_pts = 1; % Dummy
end
util.disp ( ' ' )

% Loop over all coupled states
for m=1:hamilt.coupling.n_eqs
    if hamilt.coupling.n_eqs>1
        util.disp (['State: ' hamilt.coupling.labels{m}])
    end

    % Construct filenames
    filename_re = strcat('wav_', int2str(m), '.dat');
    filename_im = strcat('wav_', int2str(m), '_im.dat');
    
    % Initialize wavefunction with zeros
    psi.dvr.grid_ND{m} = zeros(size(space.dvr.grid_ND{1}));
    
    % Check availability of input data file (real part!)
	fid = fopen(filename_re);
    if fid == -1
        util.disp (['Input data file not available : ' filename_re])
        util.disp ('Assuming real part of wave function to be zero');
        util.disp (' ')
    else
		fclose(fid);
        psi.dvr.grid_ND{m} = util.interp (...
            space.dvr.grid_ND, ...
            filename_re, ...
            psi.corr.n_pts, ...
            psi.corr.method, ...
            psi.corr.conv, ...
            1);
    end 
    
    % Check availability of input data file (imaginary part!)
	fid = fopen(filename_im);
    if fid == -1
        util.disp (['Input data file not available : ' filename_im])
        util.disp ('Assuming imaginary part of wave function to be zero');
        util.disp (' ')
    else
		fclose(fid);
        psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} + ...
            1i * util.interp (...
            space.dvr.grid_ND, ...
            filename_im, ...
            psi.corr.n_pts, ...
            psi.corr.method, ...
            psi.corr.conv, ...
            1);
    end     
    
    util.disp(' ')
end
