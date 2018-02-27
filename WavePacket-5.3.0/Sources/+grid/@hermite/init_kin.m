% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function obj = init_kin(obj, fraction, output)

global time hamilt


if nargin < 3
    output = true;
end


% Builds the matrices for the kinetic operator and propagator

%% First, let us build the matrix for the momentum operator in FBR

% The matrix form is P_{n,n+1} = conj(P_{n+1,n}) = -i sqrt((n+1)/2)
% for scaled coordinates. Unscaling induces a norm change of the
% polynomials of (m*omega)^1/4, which we have to account for here

subdiag = 1i * sqrt( (1:obj.n_pts-1) / 2 );

obj.dvr_mom = zeros(obj.n_pts) + diag(subdiag, -1) + diag(conj(subdiag), 1);
obj.dvr_mom = obj.dvr_mom * sqrt(obj.omega * obj.mass);


%% Second, transform it to DVR

% There is little use in transforming the wavefunction to
% FBR if the matrix is not diagonal there anyway. So instead of applying
% P (\Psi) = T^{-1} * M * T * \psi (T: transformation matrix, M momentum)
%, we just redefine M: M' = T^-1 * M * T
obj.dvr_mom = obj.trafo2dvr * obj.dvr_mom * obj.trafo2fbr;


%% Third: construct the kinetic energy operator

% This works independently of the representation, as they are just
% powers of the momentum operator, and the trafo DVR <=> FBR is unitary.
obj.kin = obj.dvr_mom * obj.dvr_mom / (2*obj.mass);


%% Fourth: Truncate the kinetic energy operator

% This is really tricky, because we do not have a diagonal representation
% of the operator, and just cutting off-diagonal matrix elements distorts
% the _shape_ of the operator.
%
% As a solution, we therefore have to diagonalize it, truncate the
% eigenvalues, and reconstruct the operator. The call to real() removes
% spurious complex numbers O(1e-14) that are there due to numerical
% inaccuracies.
[U, D] = eig(obj.kin);
D(D > hamilt.truncate.delta) = hamilt.truncate.delta;
obj.kin = real(U * D * inv(U));
obj.kin_max = max(D(:));


%% Fifth: create the propagator. (After long debugging...) the diag
%% functions are crucial.
obj.kin_expo = U * diag(exp(-1i * time.sub.delta * fraction * diag(D))) * inv(U);


%% Informational output

if obj.nokin
    return
end

if output
    util.disp (' ')
    util.disp ('*************************************************')
    util.disp ( ['Kinetic energy for degree of freedom: ' int2str(obj.dof)] )
    util.disp (' ')
    util.disp ('          1  (  d  ) 2                           ')
    util.disp (' T  =  - --- ( --- )                             ')
    util.disp ('         2 M ( d R )                             ')
    util.disp ('                                                 ')
    util.disp ('*************************************************')
    util.disp ( ['Mass: ' num2str(obj.mass)]       )
    util.disp (' ')
end
