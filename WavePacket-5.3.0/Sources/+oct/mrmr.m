% This function delivers a balancing transform by means of the
% Minima Realization and Model Reduction (MRMR) algorithm 
%
% J. Hahn and T. Edgar 
% Balancing Approach to Minimal Realization and 
% Model Reduction of Stable Nonlinear Systems
% Industrial & Engineering Chemistry Research 41 (9), 2204-2212 (2002)
% doi:10.1021/ie0106175 
%
% Balanced realization a for nonlinear system when the controllability
% and observability Gramians (or covariance matrices) are known.
% S: transformation matrix
% T: the inverse of the transformation matrix
% P: controllability Gramian/covariance matrix W_C
% Q: observability Gramian/covariance matrix W_O
%
% Adapted from web page of the Hahn research group at Rensselaer
% http://homepages.rpi.edu/~hahnj/Model_Reduction/index.html
% Copyright (C) 2006-2012 Chuili Sun & Juergen Hahn 
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License 

function [T, S, Sigma] = mrmr (P,Q) % Note the exchange of S and T

util.disp ('Minima Realization and Model Reduction (MRMR)')
util.disp ('=============================================')
util.disp ('   ')

n = size(P,1);
n_P = rank(P);

offdiag = zeros(n,n);
for i = 1:n
   offdiag(i,n-i+1) = 1;
end

% Schur decomposition of controllability Gramian/covariance , see Eq. (21).
% Returns unitary matrix U and upper triangular matrix T so that P = U*T*U'
% Since T is triangular, its diagonal entries contain the eigenvalues of P.
[U,T] = schur(P); 
if T(1,1)<T(n,n)
   U = U*offdiag;
end
dia = diag(T);
if dia(1)<dia(n)
   dia = offdiag*dia;
end
invdia = [sqrt(inv(diag(dia(1:n_P)))) zeros(n_P,n-n_P);zeros(n-n_P,n_P) eye(n-n_P,n-n_P)];
V = U';
T1 = invdia*V;
% T1*P*T1';

% Apply previous transformation also to observability Gramian, see Eq. (22)
Qtrans = T1'\Q/T1;
Q11 = Qtrans(1:n_P,1:n_P);

% Schur decomposition of (transformed!) observability, see Eq. (23)
[U1,~] = schur(Q11);
U1 = U1';
sigmasquared = U1*Q11*U1';

% Second transformation, see Eq. (24)
n_S = rank(sigmasquared);
T2 = inv([U1 zeros(n_P,n-n_P);zeros(n-n_P,n_P) eye(n-n_P,n-n_P)])';

% Apply transformation with T1, T2 to original observability, see Eq. (25)
Qtrans = inv(T2')*inv(T1')*Q*inv(T1)*inv(T2);
Q121 = Qtrans(1:n_S,n_P+1:n);

% Third transformation, see Eq. (26)
T3 = inv([eye(n_P) zeros(n_P,n-n_P); ...
    -Q121'/sigmasquared(1:n_S,1:n_S) zeros(n-n_P,n_P-n_S) eye(n-n_P)])';

% Apply transformation with T1, T2, T3 to original observability, see (27)
Qtrans = inv(T3')*inv(T2')*inv(T1')*Q*inv(T1)*inv(T2)*inv(T3);

% Schur decopmposition for last columns/rows of transformed system
Qt = Qtrans(n_P+1:n,n_P+1:n);
[U2,~] = schur(Qt);
U2 = U2';

% Fourth transformation, see Eq. (29)
T4 = inv([sigmasquared(1:n_S,1:n_S)^-0.25 zeros(n_S,n-n_S); ...
    zeros(n_P-n_S,n_S) eye(n_P-n_S) zeros(n_P-n_S,n-n_P); ...
    zeros(n-n_P,n_P) U2])';

T = T4*T3*T2*T1;                         % transformation matrix (original: T)
S = inv(T);                              % inverse of transformation matrix

% Transformed controllability with Sigma_1, I, 0, 0 on diagonal, see (31)
% Ptrans = T*P*T';                             % balanced controllability Gramian

% Transformed observability with Sigma_1, 0, Sigma_3, 0 on diagonal,  see (32)
% Qtrans = S'*Q*S;                             % balanced observability Gramian

% Singular values:
% Sigma=sqrt(abs(diag(Ptrans*Qtrans))); 
Sigma=sqrt(abs(diag(T*P*Q*S)));

% Sigma_1: controllable and observable
util.disp ([int2str(n_S) ' states are controllable and observable'])
% Sigma_2: controllable but not observable
util.disp ([int2str(n_P-n_S) ' states are controllable but not observable'])
% Sigma_3: observable but not controllable
% Sigma_4: neither observable nor controllable
util.disp('   ')


