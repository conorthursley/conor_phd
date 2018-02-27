% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2009 Ulf Lorenz
%
% see the README file for license details.

function NaI
global space hamilt

util.disp('  ')
util.disp('***************************************************')
util.disp('Dipole moments of the lowest two diabatic states of')
util.disp('NaI. The explicit formulas are given in            ')
util.disp('J.Chem.Soc.,Faraday Trans. 93:977.                 ')
util.disp('                                                   ')
util.disp('Values for the overlap integrals are tabulated in  ')
util.disp('J.Chem.Phys. 17:1248                               ')
util.disp('***************************************************')


mu_Na = 0.7333;         % Slater exponent for Natrium 3s orbital
mu_I  = 1.9;            % Slater exponent for Iodine 5p orbital

% If there is only one degree of freedom, we consider it as the
% internuclear distance R, for two DOF, we behave the same as
% dip.spherical_2D, i.e. fill only hamilt.d_x.grid_ND for the
% parallel dipole moment. In the latter case, we also require a
% variable R_dof set that tells us which degree of freedom is the
% radial one.

if space.size.n_dim > 2
    util.error('NaI dipole moment can only handle 1-2 dimensions')
end

if hamilt.coupling.n_eqs ~= 2
    util.error('NaI dipole moment is only for 2 coupled states')
end

if space.size.n_dim == 1
	hamilt.dip.R_dof = 1;
end

if space.size.n_dim == 2
	if ~isfield(hamilt.dip, 'R_dof') || isempty(hamilt.dip.R_dof)
		util.error('NaI dipole moment requires degree of freedom to be set.')
	end

	L_dof = 1;
	if hamilt.dip.R_dof == 1
		L_dof = 2;
	end
	angleval = space.dvr.grid_ND{L_dof};
end

Rval = space.dvr.grid_ND{hamilt.dip.R_dof};

[p, t] = convert(Rval, mu_Na, mu_I);

hamilt.d_x.grid_ND{1,1} = Rval;
hamilt.d_x.grid_ND{2,2} = Rval/2 .* mod_overlap(p,t) ./ (1 + overlap(p,t).^2);
hamilt.d_x.grid_ND{1,2} = Rval .* (overlap(p,t) + mod_overlap(p,t))...
            ./ sqrt(2 * (1+overlap(p,t).^2));

if space.size.n_dim == 2
    hamilt.d_x.grid_ND{1,1} = hamilt.d_x.grid_ND{1,1} .* angleval;
    hamilt.d_x.grid_ND{2,2} = hamilt.d_x.grid_ND{2,2} .* angleval;
    hamilt.d_x.grid_ND{1,2} = hamilt.d_x.grid_ND{1,2} .* angleval;
end

% No dipole moments along y
hamilt.d_y.grid_ND{1,1} = [];
hamilt.d_y.grid_ND{2,2} = [];
hamilt.d_y.grid_ND{1,2} = [];



function [pval, tval] = convert(R, mu_a, mu_b)
% Converts the (array of) internuclear distance, and the two slater
% exponents into the dimensionless parameters for calculating the
% overlap integral
pval = 0.5 * (mu_a + mu_b) * R;
tval = (mu_a-mu_b) / (mu_a + mu_b);




function ol = mod_overlap(p, t)
% Calculates the dipole matrix element between a STO 3s and 5p\sigma
% orbital. It occurs in equation 3 in the appendix of
% J.Chem.Soc.,Faraday Trans. 93:977, and can be evaluated using 
% J.Chem.Phys 17:1248. it turns out that it is the same as overlap(),
% but with each A/B subscript incremented by 1.
pot = p.*t;
ol = -A(8,p).*B(2, pot) + A(7,p).*B(1,pot) + 3*A(6,p).*B(4,pot) ...
     - 3*A(5,p).*B(3,pot) - 3*A(4,p).*B(6,pot) + 3*A(3,p).*B(5,pot) ...
     + A(2,p).*B(8,pot) - A(1,p).*B(7,pot);
ol = ol .* 1/(960*sqrt(42)) .* p.^8 .* (1+t).^(7/2) .* (1-t).^(9/2);



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
