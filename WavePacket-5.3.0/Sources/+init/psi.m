%------------------------------------------------------------------------------
%
% Set up the initial wave function on the different electronic states.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function psi 
global hamilt psi space

if ~isfield(psi,'init')
    psi.init = [];
end

% Null the wave function first.
for m = 1:hamilt.coupling.n_eqs
	psi.dvr.grid_ND{m} = zeros(size(space.dvr.weight_ND));
end


%% Get the wave functions
if isfield(psi, 'corr') && isfield(psi.corr, 'handle')
    % create a fully correlated wave function
    feval(psi.corr.handle);
else
    % create the wave function as direct poduct of 1D functions
    util.disp (' ')
    util.disp ('*************************************************')
    util.disp ('Initial wave function as direct product          ')
    util.disp (' ')
    util.disp ('             N                                   ')
    util.disp (' Psi(R) =  Prod   Psi  (R )                      ')
    util.disp ('            i=1      i   i                       ')
    util.disp (' ')
    util.disp ('where the product runs over all degrees of       ')
    util.disp ('freedom.                                         ')
    util.disp ('*************************************************')

    psi.dvr.grid_ND{1} = ones(size(space.dvr.grid_ND{1}));

    for k = 1:space.size.n_dim
        psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} .* feval ( psi.dof{k}.handle, k );
    end
end


%% Wavefunctions initially populated according to individual coefficients
if hamilt.coupling.n_eqs>1 && isfield(psi.init, 'coeffs') && ~isempty(psi.init.coeffs)

    % Check number of coefficients
    if length(psi.init.coeffs) ~= hamilt.coupling.n_eqs
        util.error('Wrong number of initial coefficients')
    end

    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ('Initial coefficients/populations')
    util.disp ('*******************************************************')
    util.disp ( [ '(A)diabatic representation    : ' psi.init.representation     ] )
    util.disp ( [ 'Coefficients of wavefunctions : ' num2str(psi.init.coeffs   ) ] )
    util.disp ( [ 'Populations of wavefunctions  : ' num2str(psi.init.coeffs.^2) ] )
    
    % Diabatic initial states
    for m=hamilt.coupling.n_eqs:-1:1
        psi.dvr.grid_ND {m} = psi.dvr.grid_ND{1} * psi.init.coeffs(m);
    end

    % If desired, transform initial state from the adiabatic representation
    if strcmpi(psi.init.representation,'adi')
        ket.diabatic;
    end
elseif hamilt.coupling.n_eqs > 1
    util.disp(' ')
    util.disp('*******************************************************')
    util.disp('Initial populations taken as set up by init function.  ')
    util.disp('*******************************************************')

    % has to be separated, otherwise isempty(...) gives an error
    for m = 1:min(hamilt.coupling.n_eqs, numel(psi.dvr.grid_ND))
        if isempty(psi.dvr.grid_ND{m})
            psi.dvr.grid_ND{m} = zeros(size(space.dvr.grid_ND{1}));
        end
    end
    for m = numel(psi.dvr.grid_ND)+1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = zeros(size(space.dvr.grid_ND{1}));
    end
end

%% Normalize the wavefunction
if ~isfield(psi.init, 'norm') || psi.init.norm == true
    norm2 = 0;
    for m = 1:hamilt.coupling.n_eqs
        norm2 = norm2 + sum ( abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:) );
    end

    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} / sqrt ( norm2 );
    end
end

%% Save initial wavefunction for later use (->autocorrelation)
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.init.ND{m} = psi.dvr.grid_ND{m};
end
