%------------------------------------------------------------------------------
%
% creates a transformation matrix U which reshapes a vector X 
% built from a columnwise vectorized matrix A (i.e. X=A(:)) 
% into a vector which stems from the same matrix A, but the sequence 
% of the elements is chosen in the following way:
% 1st: main diagonal of A
% 2nd: first diagonal above main diagonal
% 3rd: first diagonal below main diagonal
% 4th: second diagonal above main diagonal
% 5th: second diagonal below main diagonal
% etc.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function U =cw2df(dim)

% A transformation such as the one sought here can be easily and transparently
% constructed in three steps.

%% 1. We create the matrix representation of our density operator, where each
%     element contains the index in the old (columnwise) representation.
indexes = reshape( 1:dim^2, dim, dim );


%% 2. We create a new vector that contains in each element the index of the
%     corresponding element in the old order. First, the diagonals, then all
%     the side diagonals.
mapping = zeros(dim^2, 1);
mapping(1:dim) = diag(indexes);

start = dim+1;
for k = 1:dim-1
	num = dim - k;

	mapping(start:start+num-1) = diag(indexes, k);
	mapping(start+num:start+2*num-1) = diag(indexes, -k);

	start = start + 2*num;
end


%% 3. From this vector, we can now easily construct the transformation matrix
U = zeros(dim^2);

for index = 1:length(U)
	U(index, mapping(index)) = 1;
end
