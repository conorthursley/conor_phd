%------------------------------------------------------------------------
%
% This function creates the initial state as a bound eigenstate of a Morse
% oscillator.
% Not using the standard recursion to generate Laguerre polynomials as im-
% plemented in math.laguerre but using instead the modified version from 
% J. P. Dahl, M. Springborg,  J. Chem. Phys. 88, 4535-47 (1988)
%
% Potential function
%
%      V(r) = d_e*(1-exp(-alf*(r-r_e)))^2
%
% m_r     - (reduced) mass
% d_e     - dissociation energy
% r_e     - equilibrium position
% alf     - range parameter (steepness)
% n_q     - quantum number of the desired eigenfunction
%
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Burkhard Schmidt
%               2009 Ulf Lorenz
%
% see the README file for license details.

function init_grid = morse(dir)

global space psi

util.disp (' ')
util.disp ('*******************************************************')
util.disp ( ['Initial wavefunction for DOF :' int2str(dir)] )
util.disp ('                                                       ' )
util.disp('Morse oscillator eigenstate')

% Set default values. This simplifies live further down
if ~isfield(psi.dof{dir}, 'm_r')
        psi.dof{dir}.m_r = space.dof{dir}.mass;
end
if ~isfield(psi.dof{dir}, 'r_e')
        psi.dof{dir}.r_e = 0;
end

% Check input values
if (psi.dof{dir}.d_e<=0)
    util.error ('Dissociation energy must be positive')
end

if (psi.dof{dir}.r_e<0)
    util.error ('Equilibrium distance must be positive')
end

if (psi.dof{dir}.alf<=0)
    util.error ('Range parameter must be positive')
end
    
if (psi.dof{dir}.m_r<=0)
    util.error ('(Reduced) mass must not be positive')
end
    
% Harmonic frequency and anharmonicity constant
omega  = psi.dof{dir}.alf ...
    * sqrt ( 2*psi.dof{dir}.d_e/psi.dof{dir}.m_r );
chi    = omega^2 / ( 4*psi.dof{dir}.d_e );
lambda = omega / (2*chi); % = sqrt(2*d_e*m_r)/a_e

% Number of bound states (counting from zero!)
n_b = floor (lambda - 0.5);

% Output
util.disp ('  ')
util.disp (['Dissociation energy     : ' num2str(psi.dof{dir}.d_e)])
util.disp (['Equilibrium position    : ' num2str(psi.dof{dir}.r_e)])
util.disp (['Range parameter (alfa)  : ' num2str(psi.dof{dir}.alf)])
util.disp (['(Reduced) mass          : ' num2str(psi.dof{dir}.m_r)])
util.disp ('  ')
util.disp (['Number of bound states  : ' int2str(n_b+1)])
util.disp (['Harmonic frequency      : ' num2str(omega)])
util.disp (['Anharmonicity           : ' num2str(chi)])


% Check quantum number
if ~isfield(psi.dof{dir},'n_q')
    psi.dof{dir}.n_q = 0;
end

if (psi.dof{dir}.n_q<0)
    util.error ('quantum number must not be negative')
end

if (psi.dof{dir}.n_q>n_b)
    util.error ('quantum number exceeds number of bound states')
end

% Bound state energy
e_n = omega*(psi.dof{dir}.n_q+0.5) ...
    - chi*(psi.dof{dir}.n_q+0.5)^2 ...
    - psi.dof{dir}.d_e;
util.disp ('   ')
util.disp (['Quantum number (eigenstate) : ' int2str(psi.dof{dir}.n_q)])
util.disp (['Energy                      : ' num2str(e_n)])


% Transformed arguments; parameters
y = psi.dof{dir}.alf * ( space.dvr.grid_ND{dir}-psi.dof{dir}.r_e ); % Eq. 33
xi = 2*lambda*exp(-y); % Eq. 39
s = 2*lambda - 2*psi.dof{dir}.n_q - 1; % Eq. 46

% Recursion 
psi.n = ones(size(space.dvr.grid_ND{dir}));
if (psi.dof{dir}.n_q>0)
    psi.n1=zeros(size(space.dvr.grid_ND{dir}));
    for n=1:psi.dof{dir}.n_q
        psi.n2=psi.n1;
        psi.n1=psi.n;
        psi.n = ((2*n+s-1-xi).*psi.n1 - sqrt((n-1)*(n+s-1))*psi.n2); % Eq. 49
        psi.n = psi.n / sqrt(n*(s+n));
    end
end

% Normalization from Dahl's paper works only for nq=0 ?!?
% psi.n = psi.n .* xi.^(s/2) .* exp(-xi/2) .* sqrt(alf*s*gamma(n_q+1)/gamma(2*lambda-n_q));

% Tricky normalization (thanks to Marius Lewerenz, Marne le Vallee) 
scalog = ( log(xi)*s - xi + log(psi.dof{dir}.alf*s) - gammaln(s+1) ) / 2 ;
init_grid = psi.n .* exp(scalog); 
