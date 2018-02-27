% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%
% see the README file for license details.

function NaI
global space hamilt

util.disp('  ')
util.disp('***************************************************')
util.disp('Potential energy surfaces of the lowest two states ')
util.disp('of NaI. The explicit formulas are given in         ')
util.disp('J.Chem.Soc.,Faraday Trans. 93:977,                 ')
util.disp('values for the overlap integrals are tabulated in  ')
util.disp('J.Chem.Phys. 17:1248                               ')
util.disp('***************************************************')

%% Parameters

ion_I  = 0.384;     % ionisation energy Iodine
ion_Na = 0.1889;    % ionisation energy Natrium
aff_I  = 0.1125;    % electron affinity Iodine
aff_Na = 0.0202;    % electron affinity Natrium (J.Phys.Chem.Ref.Data 14:731)

mu_I   = 1.9;       % Slater exponent Iodine 5p STO
n_I    = 5;         % 5 orbital
mu_Na  = 0.755;     % Slater exponent Na 3s STO
n_Na   = 3;         % 3 orbital

paramA = 101.426;   % parameter A of the core potential
paramB = 2.9995;    % parameter B of the core potential
paramrho = 0.6593;  % parameter rho of the core potential

CvdW   = 18.911;    % van-der-Waals coefficient
sumpol = 46.146;    % sum of the polarizabilities; normalized such that E = sum/(2*R^4)
prodpol= 119.476;   % product of the polarizabilities; normalized such that E=2prod/R^7

%% Do some consistency check. If there is one degree of freedom, it is assumed
%% to be the radial one. For two degrees of freedom, the radial degree of
%% freedom has to be supplied via a parameter.

if hamilt.coupling.n_eqs ~= 2
    util.error('NaI potentials are only for two coupled states.')
end

if space.size.n_dim > 2
    util.error('NaI potentials are only for up to two degrees of freedom.')
end

if space.size.n_dim == 1
	hamilt.pot.R_dof = 1;
end

if space.size.n_dim == 2  ...
		&& (~isfield(hamilt.pot, 'R_dof') || isempty(hamilt.pot.R_dof))
	util.error('NaI potential requires radial degree of freedom ("R_dof").')
end

%% Calculate ionic (1) and covalent potential (2) + transition

Rval = space.dvr.grid_ND{hamilt.pot.R_dof};
[pval,tval] = convert(Rval, mu_Na, mu_I);
Sval = overlap(pval,tval);
core = ion_I + ion_Na + (paramA + (paramB./Rval).^8).*exp(-Rval/paramrho) + 1./Rval;

hamilt.pot.grid_ND{1,1} = core - ion_I - aff_I - 2./Rval ...
            - CvdW./(Rval.^6) - sumpol./(2*Rval.^4) - 2*prodpol./(Rval.^7);

hamilt.pot.grid_ND{2,2} = core + 1./(2*(1 + Sval.^2)) .* (...
            - 2*ion_I - 2*ion_Na - 2./Rval - Sval.^2/2 .* ...
            (3*ion_I + 3*ion_Na + aff_I + aff_Na + 2./Rval + 2*mu_I/n_I...
            + 2*mu_Na/n_Na));

hamilt.pot.grid_ND{1,2} = Sval./sqrt(2*(1+Sval.^2)) .* (2*core...
            - 2*ion_I - ion_Na - aff_I - 2./Rval - (mu_I/n_I + mu_Na/n_Na)/2);




function [pval, tval] = convert(R, mu_a, mu_b)
% Converts the (array of) internuclear distance, and the two slater
% exponents into the dimensionless parameters for calculating the
% overlap integral
pval = 0.5 * (mu_a + mu_b) * R;
tval = (mu_a-mu_b) / (mu_a + mu_b);




function ol = overlap(p, t)
% Calculates the overlap between a STO 3s and 5p\sigma orbital
% depending on the dimensionless parameters p and t. The equation is
% taken from J.Chem.Phys. 17:1248
pot = p .* t;
ol = -A(7,p).*B(1, pot) + A(6,p).*B(0,pot) + 3*A(5,p).*B(3,pot) ...
     - 3*A(4,p).*B(2,pot) - 3*A(3,p).*B(5,pot) + 3*A(2,p).*B(4,pot) ...
     + A(1,p).*B(7,pot) - A(0,p).*B(6,pot);
ol = ol .* 1/(960*sqrt(42)) .* p.^8 .* (1+t).^(7/2) .* (1-t).^(9/2);



function retval = B(k, pt)
% Calculates the function B_k as in J.Chem.Phys. 17:1248
% B_k(p*t) = k! * (exp(-pt) * sum_{n=1}^{k+1} 1/((pt)^n * (k-n+1)!)
%            - exp(pt) * \sum_{n=1}^{k+1} (-1)^{k-n} / ((pt)^n * (k-n+1)!)
sum1 = 0;
sum2 = 0;

for n = 1:k+1
    fac = factorial(k-n+1);
    sum1 = sum1 + 1./(pt.^n * fac);
    sum2 = sum2 + (-1)^(k-n) ./ (pt.^n * fac);
end

retval = - factorial(k) * (exp(-pt).*sum1 + exp(pt).*sum2);




function retval = A(k, p)
% Calculates the function A_k(p) as in J.Chem.Phys. 17:1248
% A_k(p) = exp(-p) * k! * sum_{n=1}^{k+1} 1/(p^n * (k-n+1)!)
retval = 0;

for n = 1:k+1
    retval = retval + 1./(p.^(n) * factorial(k-n+1));
end

retval = retval .* exp(-p) * factorial(k);
