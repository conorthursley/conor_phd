%-------------------------------------------------------------------------
%
% Operator splitting
% ==================
%
% Propagate wavefunction (psi) by time.sub.n substeps of size time.sub.delta
% by (repeated action of) Trotter or Strang splitting using FFT methods
% for a Hamiltonian operator of the type 
%
%                                        F^2(t)          (       d)
% H = V (R) - F(t) * ( D (R) + D (R) ) - ------ P(R) + T (R, -i --)
%                       p       t          2             (      dR)
%
%
% with kinetic operator T and potential V (R) energy (may be matrix-valued) 
% and where the interaction between an (optional) external electric 
% field F(t) and the dipole moment D(R) and/or the polarizability P(R)
% of the system is treated semiclassically. Here, the dipole matrix D 
% has been split into its diagonal (D_p) and offdiagonal (D_t) parts
% containing permanent and transition dipoles, respectively. In most
% cases, only one of the two is expected to play a role, depending on 
% the light frequencies under consideration. Moreover, note that F(t) 
% as well D_p(R), D_t(R), and P(R) can have two cartesian components 
% (along x,y) corresponding to different polarization directions.
%
% (1) Exact (for time-independent Hamiltonian)
%
% psi(t+tau) = exp(-i*H*tau  ) * psi(t) 
%
% (2) Trotter splitting
%
% psi(t+tau) = exp( -i*         V *tau ) 
%            * exp( +i*F(t)    *Dp*tau ) 
%            * exp( +i*F(t)    *Dt*tau ) 
%            * exp( +i*F^2(t)/2*P *tau ) 
%            * exp( -i*         T *tau )
%            * psi(t)                    + O(tau^2)
%
% (3) Strang splitting
%
% psi(t+tau) = exp( -i*         V *tau/2 ) 
%            * exp( +i*F(t)    *Dp*tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( -i*         T *tau   ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 )
%            * exp( +i*F(t)    *Dp*tau/2 )
%            * exp( -i*         V *tau/2 ) 
%            * psi(t)                    + O(tau^3)
%
%
% Implementation for a single Schrödinger equation
% ------------------------------------------------
%
% The above Trotter and Strang formulae are straightforward to evaluate for
% FFT grids (plane wave DVR):
%   Operator V, and hence exp(V), are diagonal in position representation
%   Operator D, and hence exp(D), are diagonal in position representation
%   Operator T, and hence exp(T), are diagonal in momentum representation
% Hence, two FFTs have to be performed for every timestep to transform
% from position to momentum representation and back. Similar arguments
% apply for other DVR methods.
%
% Implementation for coupled Schrödinger equation
% -----------------------------------------------
%
% In case of coupled Schrödinger equations, the potential energy V as well
% as the dipole moment operator D may be represented by (real, symmetric) 
% matrices. Hence, for each spatial discretization point, the exponential 
% of these matrices has to be calculated. This can be achieved either by
% Matlab's expm function or by the method of eigenvectors and eigenvalues, 
% where the accuracy is determined by the condition of the matrix.
% 
% For the case of two equations, the matrix exponential can be calculated 
% analytically, see  Eq. (11.204) on page 316 of David Tannor's book. 
% 
%     ( alpha   beta )
% V = (              )
%     ( beta^* gamma )
%
% delta = (alfa-gamma)/2
%   eta = (alfa+gamma)/2
%   rho = sqrt(delta^2+beta^2)
%
% expm(-iV*tau) 
%   = exp(-i*eta*tau) *
%   * [              cos(tau*rho)*I 
%       -i*delta/rho*sin(tau*rho)*S3
%       -i* beta/rho*sin(tau*rho)*S1  ]
%
% where I=(1 0; 0 1), S1 = (0 1; 1 0), S3 = (1 0; 0 -1)
%
% 
% References
% ==========
%
% Original version: 
%     J. A. Fleck et al.,                  Appl. Phys.,   10,  129 (1976)
%     M. D. Feit, J. A. Fleck, A. Steiger, J. Comp. Phys. 47,  412 (1982) 
% Coupled equations: 
%     J. Alvarellos and H. Metiu,          J. Chem. Phys. 88, 4957 (1988)
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function splitting ( step )
global psi space time hamilt

% First step only: Initialize propagations
if step==1
    
    util.disp('  ')
    util.disp('*************************')
    util.disp('Operator splitting method')
    util.disp('*************************')
    
    if ~isfield (time.propa,'order')
        time.propa.order = 3;
    end
    
    switch time.propa.order
        case (2)
            util.disp ( 'Trotter splitting' )
            fraction = 1;
        case (3)
            util.disp ( 'Strang splitting' )
            fraction = 1/2;
        otherwise
            util.error ( 'Wrong choice of splitting order' )
    end
    util.disp ( ' ' )
    
    % Propagator for kinetic energy
    for k = 1:space.size.n_dim-1
        space.dof{k} = init_kin(space.dof{k}, fraction, false);
    end
    
    % External kinetic energies; Since they are typically the most
    % time-consuming operators, we put them in the middle of the split operator
    % If none are specified, take the last grid-internal kinetic energy
    % in the middle.
    if isfield(hamilt, 'kin')
        for n = 1:length(hamilt.kin)-1
            hamilt.kin{n} = init_kin(hamilt.kin{n}, fraction, false);
        end
        hamilt.kin{end} = init_kin(hamilt.kin{end}, 1, false);
        space.dof{end} = init_kin(space.dof{end}, fraction, false);
    else
        space.dof{end} = init_kin(space.dof{end}, 1, false);
    end
   
	% Propagator for the potential energy
    pot_init (fraction)
    
    % Propagator for (permanent and/or transition) dipoles, polarizabilities
    if time.efield.n_pulse > 0
        perm_init (fraction);
        trans_init (fraction);
        pol_init (fraction);
    end

% After first step: Perform propagations
else

    % Loop over N sub-steps: 
    for k = 1:time.sub.n

        % Potential energy: Full (Trotter) or half (Strang) sub-step
        pot_propa;

        % Permanent/transition dipoles, polarizabilities: Full (Trotter) or half (Strang) sub-step
        if time.efield.n_pulse > 0
            e_x = time.efield.short.x(k);
            e_y = time.efield.short.y(k);
            perm_propa  (e_x, e_y);
            trans_propa (e_x, e_y, k==1);
            pol_propa   (e_x, e_y);
        end


        % Kinetic energy: Propagate each DOF
        for l = 1:space.size.n_dim
            kinetic_exp(space.dof{l});
        end

        % Propagate each external kinetic energy
        if isfield(hamilt, 'kin')
            for n = 1:length(hamilt.kin)
                apply_exp(hamilt.kin{n});
            end
        end

        % Strang only 
        if time.propa.order==3 
            
            % Inverse order of the external kinetic energies
            % And apply the last grid kinetic operator as well.
            if isfield(hamilt, 'kin')
                for n = length(hamilt.kin)-1:-1:1
                    apply_exp(hamilt.kin{n});
                end
                kinetic_exp(space.dof{end});
            end

            % Apply all but the last grid kinetic energy operator
            % (it has either been applied before, or it was the
            % innermost split operator).
            for l = space.size.n_dim-1:-1:1
                kinetic_exp(space.dof{l});
            end

            % Transition/permanent dipoles, polarizabilities
            if time.efield.n_pulse>0
                e_x = time.efield.short.x(k+1);
                e_y = time.efield.short.y(k+1);

                pol_propa   (e_x, e_y);
                trans_propa (e_x, e_y, true);
                perm_propa  (e_x, e_y);
            end

            % Potential energy
            pot_propa;
        
        end
        
        % Autocorrelation function
        time.acf.short(k) = 0;
        for m = 1:hamilt.coupling.n_eqs
            time.acf.short(k) = time.acf.short(k) + ...
                sum ( conj(psi.dvr.init.ND{m}(:)).*psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:) );
        end

    end

end



%------------------------------------------------------------------
% Initialize propagation associated with polarizabilities (vector)
%------------------------------------------------------------------
function pol_init (fraction)
global hamilt time 

% Detect any polarizabilities along x
hamilt.p_x.toggle = false;
if isfield (hamilt,'p_x')
    if isfield (hamilt.p_x,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.p_x.grid_ND{m,m} )
                hamilt.p_x.toggle = true;
            end
        end
    end
end

% Detect any polarizabilities along y
hamilt.p_y.toggle = false;
if isfield (hamilt,'p_y')
    if isfield (hamilt.p_y,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.p_y.grid_ND{m,m} )
                hamilt.p_y.toggle = true;
            end
        end
    end
end

% Polarizabilities along x
if hamilt.p_x.toggle
    util.disp(['Initialize x-polarizabilities propagator for time step fraction: ', num2str(fraction)])
    for m=1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.p_x.grid_ND{m,m} )
            hamilt.p_x.prod{m,m}  = -1i * time.sub.delta*fraction * hamilt.p_x.grid_ND{m,m};
        end
    end
end

% Polarizabilities along y
if hamilt.p_y.toggle
    util.disp(['Initialize y-polarizabilities propagator for time step fraction: ', num2str(fraction)])
    for m=1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.p_y.grid_ND{m,m} )
            hamilt.p_y.prod{m,m}  = -1i * time.sub.delta*fraction * hamilt.p_y.grid_ND{m,m};
        end
    end
end
 
%---------------------------------------------------------------
% Perform propagation associated with polarizabilities (vector)
%---------------------------------------------------------------
function pol_propa ( e_x, e_y )
global psi hamilt

if hamilt.p_x.toggle && abs(e_x)>0
    for m = 1:hamilt.coupling.n_eqs
        if  ~isempty ( hamilt.p_x.prod{m,m} )
            psi.dvr.grid_ND{m} = exp(hamilt.p_x.prod{m,m} * (-e_x.^2/2)) .* psi.dvr.grid_ND{m};
        end
    end
end

if hamilt.p_y.toggle && abs(e_y)>0
    for m = 1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.p_y.prod{m,m} )
            psi.dvr.grid_ND{m} = exp(hamilt.p_y.prod{m,m} * (-e_y.^2/2)) .* psi.dvr.grid_ND{m};
        end
    end
end

%------------------------------------------------------------------
% Initialize propagation associated with permanent dipoles (vector)
%------------------------------------------------------------------
function perm_init (fraction)
global hamilt time 

% Detect any permanent dipole moments along x
hamilt.d_x.perm = false;
if isfield (hamilt,'d_x')
    if isfield (hamilt.d_x,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.d_x.grid_ND{m,m} )
                hamilt.d_x.perm = true;
            end
        end
    end
end

% Detect any permanent dipole moments along y
hamilt.d_y.perm = false;
if isfield (hamilt,'d_y')
    if isfield (hamilt.d_y,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.d_y.grid_ND{m,m} )
                hamilt.d_y.perm = true;
            end
        end
    end
end

% Permanent dipole moment along x
if hamilt.d_x.perm
    util.disp(['Initialize permanent x-dipole propagator for time step fraction: ', num2str(fraction)])
    for m=1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.d_x.grid_ND{m,m} )
            hamilt.d_x.prod{m,m}  = -1i * time.sub.delta*fraction * hamilt.d_x.grid_ND{m,m};
        end
    end
end

% Permanent dipole moment along y
if hamilt.d_y.perm
    util.disp(['Initialize permanent y-dipole propagator for time step fraction: ', num2str(fraction)])
    for m=1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.d_y.grid_ND{m,m} )
            hamilt.d_y.prod{m,m}  = -1i * time.sub.delta*fraction * hamilt.d_y.grid_ND{m,m};
        end
    end
end
 
%---------------------------------------------------------------
% Perform propagation associated with permanent dipoles (vector)
%---------------------------------------------------------------
function perm_propa ( e_x, e_y )
global psi hamilt

if hamilt.d_x.perm && abs(e_x)>0
    for m = 1:hamilt.coupling.n_eqs
        if  ~isempty ( hamilt.d_x.prod{m,m} )
            psi.dvr.grid_ND{m} = exp(hamilt.d_x.prod{m,m} * (-e_x)) .* psi.dvr.grid_ND{m};
        end
    end
end

if hamilt.d_y.perm && abs(e_y)>0
    for m = 1:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.d_y.prod{m,m} )
            psi.dvr.grid_ND{m} = exp(hamilt.d_y.prod{m,m} * (-e_y)) .* psi.dvr.grid_ND{m};
        end
    end
end

%-------------------------------------------------------------------
% Initialize propagation associated with transition dipoles (matrix)
%-------------------------------------------------------------------
function trans_init (fraction)
global hamilt space time

% Detect any transition dipole moment along x
hamilt.d_x.trans = false;
if isfield (hamilt,'d_x')
    if isfield (hamilt.d_x,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            for n=m+1:hamilt.coupling.n_eqs
                if  ~isempty(hamilt.d_x.grid_ND{m,n})
                    hamilt.d_x.trans = true;
                end
            end
        end
    end
end

% Detect any transition dipole moment along y
hamilt.d_y.trans = false;
if isfield (hamilt,'d_y')
    if isfield (hamilt.d_y,'grid_ND')
        for m=1:hamilt.coupling.n_eqs
            for n=m+1:hamilt.coupling.n_eqs
                if ~isempty(hamilt.d_y.grid_ND{m,n})
                    hamilt.d_y.trans = true;
                end
            end
        end
    end
end


switch hamilt.coupling.n_eqs
 
    % No transition dipoles for one equation (nothing to be done)
    case 1
               
    % Analytical eigen/values/vectors for two coupled equations
    % [V,D]=eig(-i*tau*[0 d12;d12 0])
    case 2
 
        if hamilt.d_x.trans
            util.disp(['Initialize analytical x-transition dipole propagator for time step fraction: ', num2str(fraction)])
            hamilt.d_x.eig_vals{1} = + 1i * time.sub.delta*fraction * hamilt.d_x.grid_ND{1,2};
            hamilt.d_x.eig_vals{2} = - 1i * time.sub.delta*fraction * hamilt.d_x.grid_ND{1,2};
            hamilt.d_x.eig_vecs{1,1} = - ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_x.eig_vecs{1,2} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_x.eig_vecs{2,1} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_x.eig_vecs{2,2} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
       end
        if hamilt.d_y.trans
            util.disp(['Initialize analytical y-transition dipole propagator for time step fraction: ', num2str(fraction)])
            hamilt.d_y.eig_vals{1} = + 1i * time.sub.delta*fraction * hamilt.d_y.grid_ND{1,2};
            hamilt.d_y.eig_vals{2} = - 1i * time.sub.delta*fraction * hamilt.d_y.grid_ND{1,2};
            hamilt.d_y.eig_vecs{1,1} = - ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_y.eig_vecs{1,2} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_y.eig_vecs{2,1} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
            hamilt.d_y.eig_vecs{2,2} = + ones(size(hamilt.pot.grid_ND{1,1}))/sqrt(2);
        end
                
    % Numerical eigen/values/vectors for more than two coupled equations 
    otherwise

        % Transition dipole moments along x direction
        if hamilt.d_x.trans

            util.disp(['Initialize numerical x-transition dipole propagator for time step fraction: ', num2str(fraction)])
 
            % Preallocate
            d_x_mat = zeros(hamilt.coupling.n_eqs);
            hamilt.d_x.eig_vals = cell(hamilt.coupling.n_eqs,1);
            hamilt.d_x.eig_vecs = cell(hamilt.coupling.n_eqs  );
            for m=1:hamilt.coupling.n_eqs
                hamilt.d_x.eig_vals{m} = zeros(size(hamilt.pot.grid_ND{1,1}));
                for n=1:hamilt.coupling.n_eqs % Full matrices: not symmetric
                    hamilt.d_x.eig_vecs{m,n} = zeros(size(hamilt.pot.grid_ND{1,1}));
                end
            end

            % Loop over all grid points
            for gr=1:space.size.n_tot
                
                % Assemble matrix of available transition dipoles
                for m=1:hamilt.coupling.n_eqs
                    for n=m+1:hamilt.coupling.n_eqs
                        if ~isempty ( hamilt.d_x.grid_ND{m,n} )
                            d_x_mat(m,n) = hamilt.d_x.grid_ND{m,n}(gr);
                            d_x_mat(n,m) = hamilt.d_x.grid_ND{m,n}(gr);
                        end
                    end % for n
                end % for m

                % Numerical eigen/values/vectors, save for later use
                [V,D] = eig(d_x_mat );
                for m=1:hamilt.coupling.n_eqs
                    % Don't put the factors in the eigenvector routine! It makes
                    % the calculation slower (complex instead of real arithmetic),
                    % and can lead to problems with non-orthonormalized eigenvectors!
                    hamilt.d_x.eig_vals{m}(gr) = D(m,m) * -1i * time.sub.delta * fraction;
                    for n=1:hamilt.coupling.n_eqs % Full matrix
                        hamilt.d_x.eig_vecs{m,n}(gr) = V(m,n);
                    end % for n
                end % for m

            end % for gr
        end % if hamilt.d_x.trans
        

        % Transition dipole moments along y direction
        if hamilt.d_y.trans
            
            util.disp(['Initialize numerical y-transition dipole propagator for time step fraction: ', num2str(fraction)])

            % Preallocate
            d_y_mat = zeros(hamilt.coupling.n_eqs);
            hamilt.d_y.eig_vals = cell(hamilt.coupling.n_eqs,1);
            hamilt.d_y.eig_vecs = cell(hamilt.coupling.n_eqs  );
            for m=1:hamilt.coupling.n_eqs
                hamilt.d_y.eig_vals{m} = zeros(size(hamilt.pot.grid_ND{1,1}));
                for n=1:hamilt.coupling.n_eqs % Full matrices: not symmetric
                    hamilt.d_y.eig_vecs{m,n} = zeros(size(hamilt.pot.grid_ND{1,1}));
                end
            end
            
            % Loop over all grid points
            for gr=1:space.size.n_tot

                % Assemble matrix of available transition dipoles
                for m=1:hamilt.coupling.n_eqs
                    for n=m+1:hamilt.coupling.n_eqs
                        if ~isempty ( hamilt.d_y.grid_ND{m,n} )
                            d_y_mat(m,n) = hamilt.d_y.grid_ND{m,n}(gr);
                            d_y_mat(n,m) = hamilt.d_y.grid_ND{m,n}(gr);
                        end
                    end % for n
                end % for m

                % Numerical eigen/values/vectors, save for later use
                [V,D] = eig(d_y_mat);
                for m=1:hamilt.coupling.n_eqs
                    hamilt.d_y.eig_vals{m}(gr) = D(m,m) * -1i * time.sub.delta*fraction;
                    for n=1:hamilt.coupling.n_eqs  % Full matrix
                        hamilt.d_y.eig_vecs{m,n}(gr) = V(m,n);
                    end % for n
                end % for m

            end % for gr
        end  % if hamilt.d_y.trans

end % switch n_eqs
        
%----------------------------------------------------------------
% Perform propagation associated with transition dipoles (matrix)
%----------------------------------------------------------------
function trans_propa (e_x,e_y,recalc)
global psi hamilt time
persistent diagonal_x diagonal_y

if hamilt.d_x.trans && abs(e_x)>0

	if recalc || isempty(diagonal_x)
		
		% HACK: If we allow complex electric fields (i.e., rotating wave
		% approximation), the eigenvectors become dependant on the phase of the
		% electric field, and we need to recalculate them over and over again.
		if time.efield.complex
			phase = e_x / abs(e_x);

			hamilt.d_x.eig_vecs{1,1} = -ones(size(hamilt.pot.grid_ND{1,1})) * phase / sqrt(2);
			hamilt.d_x.eig_vecs{2,1} = ones(size(hamilt.pot.grid_ND{1,1})) / sqrt(2);
			hamilt.d_x.eig_vecs{1,2} = ones(size(hamilt.pot.grid_ND{1,1})) * phase / sqrt(2);
			hamilt.d_x.eig_vecs{2,2} = ones(size(hamilt.pot.grid_ND{1,1})) / sqrt(2);

			for m = 1:hamilt.coupling.n_eqs
				diagonal_x{m} = exp(-abs(e_x)*hamilt.d_x.eig_vals{m});
			end
		else
			% Profiling has shown that for more realistic settings (H2+ in a laser
			% field with 3 coupled surfaces) the exponentiation here easily makes up
			% 1/4-1/3 of the computing time. By calculating it only once per
			% timestep, we can thus shave off >10%.
			for m = 1:hamilt.coupling.n_eqs
				diagonal_x{m} = exp(-e_x*hamilt.d_x.eig_vals{m});
			end
		end
	end
            
    
    % Transform to adiabatic (transformation matrix D+) 
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.new_ND{m} = zeros(size(psi.dvr.grid_ND{1}));
        for n=1:hamilt.coupling.n_eqs
            psi.dvr.new_ND{m} = psi.dvr.new_ND{m} ...
                + conj(hamilt.d_x.eig_vecs{n,m}) .* psi.dvr.grid_ND{n};
        end

        % Propagate in adiabatic representation
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} .* diagonal_x{m};
    end
    
    % Transform back to diabatic (transformation matrix D) 
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = zeros(size(psi.dvr.new_ND{1}));
        for n=1:hamilt.coupling.n_eqs
            psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} ...
                + hamilt.d_x.eig_vecs{m,n} .* psi.dvr.new_ND{n};
        end
    end
    
end

if hamilt.d_y.trans && abs(e_y)>0

	if recalc || isempty(diagonal_x)
		% HACK: If we allow complex electric fields (i.e., rotating wave
		% approximation), the eigenvectors become dependant on the phase of the
		% electric field, and we need to recalculate them over and over again.
		if time.efield.complex
			phase = e_y / abs(e_y);

			hamilt.d_y.eig_vecs{1,1} = -ones(size(hamilt.pot.grid_ND{1,1})) * phase / sqrt(2);
			hamilt.d_y.eig_vecs{2,1} = ones(size(hamilt.pot.grid_ND{1,1})) / sqrt(2);
			hamilt.d_y.eig_vecs{1,2} = ones(size(hamilt.pot.grid_ND{1,1})) * phase / sqrt(2);
			hamilt.d_y.eig_vecs{2,2} = ones(size(hamilt.pot.grid_ND{1,1})) / sqrt(2);

			for m = 1:hamilt.coupling.n_eqs
				diagonal_y{m} = exp(-abs(e_y)*hamilt.d_y.eig_vals{m});
			end
		else
			% Profiling has shown that for more realistic settings (H2+ in a laser
			% field with 3 coupled surfaces) the exponentiation here easily makes up
			% 1/4-1/3 of the computing time. By calculating it only once per
			% timestep, we can thus shave off >10%.
			for m = 1:hamilt.coupling.n_eqs
				diagonal_y{m} = exp(-e_y*hamilt.d_y.eig_vals{m});
			end
		end
	end

    % Transform to adiabatic (transformation matrix D+) 
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.new_ND{m} = zeros(size(psi.dvr.grid_ND{1}));
        for n=1:hamilt.coupling.n_eqs
            psi.dvr.new_ND{m} = psi.dvr.new_ND{m} ...
                + conj(hamilt.d_y.eig_vecs{n,m}) .* psi.dvr.grid_ND{n};
        end
            
        % Propagate in adiabatic representation
        psi.dvr.new_ND{m} = psi.dvr.new_ND{m} .* diagonal_y{m};
    end
    
    % Transform back to diabatic (transformation matrix D) 
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = zeros(size(psi.dvr.new_ND{1}));
        for n=1:hamilt.coupling.n_eqs
            psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} ...
                + hamilt.d_y.eig_vecs{m,n} .* psi.dvr.new_ND{n};
        end
    end
    
end


%----------------------------------------------------------------
% Initialize propagator associated with potential energy (matrix)
%----------------------------------------------------------------
function pot_init (fraction)
global hamilt space time
    
% Detect whether any potential coupling exists
hamilt.pot.couple = false;
for m=1:hamilt.coupling.n_eqs
    for n=m+1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.pot.grid_ND{m,n})
            hamilt.pot.couple = true;
        end
    end
end

% Diabatic potential coupling is absent
if  ~hamilt.pot.couple
    util.disp(['Initialize diagonal potential propagator for time step fraction: ', num2str(fraction)])
    for m=1:hamilt.coupling.n_eqs
        hamilt.pot.expo{m,m}  = exp ( -1i * time.sub.delta*fraction * hamilt.pot.grid_ND{m,m}   );
    end

% Potential coupling exists
else

    switch hamilt.coupling.n_eqs

        % Scalar exponential for one equation (nothing to be done)
        case 1

        % Analytic matrix exponential for two equations
        % expm(-i*tau*[eta+delta beta; beta eta-delta])
        case 2
            util.disp(['Initialize analytical potential coupling propagator for time step fraction: ', num2str(fraction)])
            delta = ( hamilt.pot.grid_ND{1,1}-hamilt.pot.grid_ND{2,2} ) /2;
            eta   = ( hamilt.pot.grid_ND{1,1}+hamilt.pot.grid_ND{2,2} ) /2;
            beta  =   hamilt.pot.grid_ND{1,2};
            rho   = sqrt (delta.^2 + beta.^2);

            ex = exp( -1i * time.sub.delta*fraction * eta );
            co = cos(       time.sub.delta*fraction * rho );
            si = sin(       time.sub.delta*fraction * rho );

            hamilt.pot.expo{1,1}  = ex .* (co - 1i*delta.*si./rho);
            hamilt.pot.expo{2,2}  = ex .* (co + 1i*delta.*si./rho);
            hamilt.pot.expo{1,2}  = ex .* (   - 1i* beta.*si./rho);

            % Special care needs to be taken at conical interaction(s)
            % alfa=gamma=eta, delta=0, beta=0, rho=0
            % Hence, potential matrix is diagonal
            ci = find(rho(:)==0); % Using linear indexing
            hamilt.pot.expo{1,1}(ci)  = ex(ci);
            hamilt.pot.expo{2,2}(ci)  = ex(ci);
            hamilt.pot.expo{1,2}(ci)  = 0;
            
        % Numerical matrix exponential for more than two coupled equations
        otherwise

            util.disp(['Initialize numerical potential coupling propagator for time step fraction: ', num2str(fraction)])
 
            % Preallocate
            pot_mat = zeros(hamilt.coupling.n_eqs);
            hamilt.pot.expo = cell(hamilt.coupling.n_eqs);
            for m=1:hamilt.coupling.n_eqs
                for n=m:hamilt.coupling.n_eqs
                    hamilt.pot.expo{m,n} = zeros(size(hamilt.pot.grid_ND{1,1}));
                end
            end

            % Loop over all points of (multidimensional) grid
            for gr=1:space.size.n_tot

                % Extract diabatic potentials from cell array
                for m=1:hamilt.coupling.n_eqs
                    pot_mat(m,m) = hamilt.pot.grid_ND{m,m}(gr);

                    % Diabatic potential coupling (Hermitian)
                    for n=m+1:hamilt.coupling.n_eqs
                        if ~isempty(hamilt.pot.grid_ND{m,n})
                            pot_mat(m,n) = hamilt.pot.grid_ND{m,n}(gr);
                            pot_mat(n,m) = hamilt.pot.grid_ND{m,n}(gr);
                        end
                    end

                end

                % Either use Matlab's expm function ...
                % exp_pot = expm ( -i * time.sub.delta*fraction * pot_mat );

                % ... or use eigen/values/vectors (from expmdemo3.m)
                [V,D] = eig(pot_mat);
                exp_pot = V * diag(exp(-1i * time.sub.delta * fraction * diag(D))) / V;

                % Store potential propagator for later use
                % Real, symmetic matrices
                for m=1:hamilt.coupling.n_eqs
                    for n=m:hamilt.coupling.n_eqs
                        hamilt.pot.expo{m,n}(gr) = exp_pot(m,n);
                    end
                end

            end % loop over grid points

    end % switch n_eqs

end % if potential coupling exists


%--------------------------------------------------------------
% Perform propagation associated with potential energy (matrix)
%--------------------------------------------------------------
function pot_propa
global psi hamilt

% Dynamics along diabatic potential energy surfaces 
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.new_ND{m} = hamilt.pot.expo{m,m} .* psi.dvr.grid_ND{m};
end

% Potential coupling
if hamilt.pot.couple
    for m = 1:hamilt.coupling.n_eqs
        for n = m+1:hamilt.coupling.n_eqs
            if ~isempty ( hamilt.pot.expo{m,n} )
                psi.dvr.new_ND{m} = psi.dvr.new_ND{m} + hamilt.pot.expo{m,n} .* psi.dvr.grid_ND{n};
                psi.dvr.new_ND{n} = psi.dvr.new_ND{n} + hamilt.pot.expo{m,n} .* psi.dvr.grid_ND{m};
            end
        end
    end
end

% Save 
psi.dvr.grid_ND = psi.dvr.new_ND;


