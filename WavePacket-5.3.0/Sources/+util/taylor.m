%*******************************************************
%
% Taylor series expansion (diagonal in N dimensions)
%
%         inf   N   F_jk           j             
% f (R) = Sum  Sum  ---- ( R  - S )  + C        
%         j=1  k=1   j!     k    k            
%
%*******************************************************

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 Burkhard Schmidt
%
% see the README file for license details.

function fr = taylor (R,S,C,F)

% Confirm input variables
util.disp ( [ '==> Reference position S      : ' num2str(S) ] )
util.disp ( [ '==> Constant offset C         : ' num2str(C) ] )

% Use constant offset to initialize function
fr = C * ones(size(R{1}));
        
% Size of coefficient matrix
[nrow,ncol]=size(F);
if ncol~=length(S)
    util.error ('==> Wrong length of Taylor coefficient vectors')
else
    util.disp (['==> Taylor expansion order    : ' int2str(nrow)])
end

% Loop over expansion orders
for j=1:nrow
    fj = factorial(j);
    
    % Tabulate expansion coefficients
    util.disp ( [ '==> ' int2str(j) '-th order coefficient(s) : ' num2str(F(j,:)) ] )
    
    % Summing up contributions from each component of position vector
    for k = 1:ncol
        fr = fr + (R{k}-S(k)).^j * F(j,k) / fj;
    end
end