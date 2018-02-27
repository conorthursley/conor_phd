%----------------------------------------------------------------------
%
% This function sets up the grids for the coordinates in position and
% "momentum" space (DVR and FBR, respectively) as well as the weights
% needed for integration in DVR.
%
%----------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function grid
global space 

% Dimensionality of problem and total number of grid points
space.size.n_dim = length(space.dof);
space.size.n_tot = 1;
for k = 1:space.size.n_dim
    space.size.n_tot = space.size.n_tot * space.dof{k}.n_pts;
end

% Some logging header
util.disp ( '**************************************************' )         
util.disp ( 'Direct product of one-dimensional grids           ' )
util.disp ( [ 'Dimensionality of problem          : ' int2str(space.size.n_dim) ] )
util.disp ( [ 'Total number of grid points        : ' int2str(space.size.n_tot) ] )
util.disp ( '**************************************************' )
util.disp ( ' ' )

% Build one-dimensional grids in DVR and FBR space as cell arrays, where
% the cell index denotes the DOF. Generates also console / log file output.
for k = 1:space.size.n_dim
    space.dof{k}.dof = k; % Tell each DVR which degree of freedom it represents.
    [space.dof{k} space.dvr.weight_1D{k} space.dvr.grid_1D{k} space.fbr.grid_1D{k}] = init_grid(space.dof{k});
    disp (space.dof{k})
end

% For multidimensional calculations, construct multi-dimensional grids 
if space.size.n_dim == 1 % Note that ndgrid(x) has the same effect as ndgrid(x,x)
    space.dvr.grid_ND{1} = space.dvr.grid_1D{1};
    space.fbr.grid_ND{1} = space.fbr.grid_1D{1};
    tmpweights{1} = space.dvr.weight_1D{1};
else 
    [space.dvr.grid_ND{1:space.size.n_dim}] = ndgrid(space.dvr.grid_1D{1:space.size.n_dim});
    [space.fbr.grid_ND{1:space.size.n_dim}] = ndgrid(space.fbr.grid_1D{1:space.size.n_dim});
    [tmpweights{1:space.size.n_dim}] = ndgrid(space.dvr.weight_1D{1:space.size.n_dim});
end

% Finally, put all the weights together to get a matrix that has to be applied
% when summing over the position grid.
if space.size.n_dim == 1
    space.dvr.weight_ND = tmpweights{1};
else
    space.dvr.weight_ND = tmpweights{1};
    for k = 2:space.size.n_dim
        space.dvr.weight_ND = space.dvr.weight_ND .* tmpweights{k};
    end
end
