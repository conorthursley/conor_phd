%--------------------------------------------------------------------------
%
% Provide input for bilinear control problems
%
% Matrix   *A*  is created from energies (and coupling to bath for LvNE)
% Vectors  *B* are created from transition dipole moments
% Matrices *N* are created from transition dipole moments
% Vectors  *C* are created from (linear) observables (LvNE)
% Matrices *D* are created from (quadratic) observables (TDSE)
% Vector *x.initial* is the initial state
% Vector *y.initial* is the corresponding output
% Vector *x.equilib* is the equilibrium state (fix point of A)
% Vector *y.equilib* is the corresponding output
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014-17 Burkhard Schmidt
%               2012 Boris Schaefer-Bung, Burkhard Schmidt, 
%                    Ulf Lorenz, Jeremy Rodriguez
%
% see the README file for license details.
%

function qm_abncd(eom)

global control bilinear 

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Load the energies, dipole, and (system-bath) coupling matrices
load ('tise');
dim=size(tise.ham, 1);

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' Provides matrices A, N, B and C, D and vectors x_i, x_e     ')
util.disp (' for use in a bilinear control problem                       ')
util.disp (' https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_abncd/')
util.disp ('                                                             ')
util.disp ('          d                                                  ')
util.disp (' control: -- x(t) = ( A + iu(t)N ) x(t) + iu(t)B             ')
util.disp ('          dt                                                 ')
util.disp ('                           T                                 ')
util.disp (' observe: y(t) = C x(t) + x (t) D x(t)                       ')
util.disp ('                                                             ')
util.disp (' see B. Schaefer-Bung, C. Hartmann, B. Schmidt, Ch. Schuette ')
util.disp (' J. Chem. Phys. 135, 014112-1-13 (2011)                      ')

switch lower(eom)
    case 'lvne'  
        % Map matrices onto vectors using the so-called tetradic (Liouville) 
        % convention introduced by Mukamel's group
        %
        % We use columnwise ordering of the matrix elements, e.g. for
        % the density rho = (rho_00, rho_10, ..., rho_01, rho_11, ...)
        % Note that this is the standard ordering of Fortran/Matlab anyway
        
        util.disp (' ')
        util.disp (' coming from the quantum Liouville-von Neumann equation')
        util.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
        util.disp (' for quantum systems interacting with electrical fields')
        util.disp ('                                                       ')
        util.disp ('  d             i                                      ')
        util.disp (' -- rho(t) = - ---- [H - F(t) mu, rho(t)] +L [rho(t)]  ')
        util.disp (' dt            hbar   0                    D           ')
        util.disp ('                                                       ')
        util.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
        util.disp (' mu is the Hermitian matrix with dipole moments        ')
        util.disp (' and F(t) is the electric field.                       ')
        util.disp ('                                                       ')
        util.disp (' with L [rho] = ... Lindblad dissipation/dephasing ... ')
        util.disp ('       D                                               ')
        util.disp ('                                                       ')
        util.disp (' <C>(t) = tr( C rho(t) )                               ')
        util.disp ('                                                       ')
        util.disp ('-------------------------------------------------------')
        util.disp (' ')
        util.disp ('Parameters of dissipative Liouvillian: Lindblad form')
        
        % Initial and equilibrium state vectors/densities for propagation
        oct.lvne (real(tise.ham))

        % set temperature; generate title string for plots
        if ~isfield(control.lvne,'temperature')
            control.lvne.temperature = 0;
        end
        util.disp (['temperature * k_B = ' num2str(control.lvne.temperature)]) 
        util.disp (' ') 
        kBT = control.lvne.temperature;
        bilinear.title = ['LvNE: \Theta =' num2str(control.lvne.temperature) ', '];
  
        % Construct omega-matrix (Bohr frequencies)
        w=zeros(dim);
        for k=1:dim
            for l=1:dim
                w(k,l)=tise.ham(k)-tise.ham(l);
                w = real(w); % to be done: delete?
            end
        end
        
        % Set up matrix A; coherent evolution from omega-matrix
        bilinear.A = diag( reshape(-1i*w, dim^2, 1) );

        % Lindblad operators for population relaxation
        if ~isfield(control,'relax')
            control.relax = [];
        end
        if ~isfield(control.relax,'model')
            control.relax.model = 'none';
        end
        
        % Construct Gamma matrix
        % convention used throughout: Gam(k,l) = Gamma_{k <- l}
        Gam = zeros(dim);
        
        switch lower(control.relax.model)
            
            % Relaxation rates from Fermi's golden rule, see Andrianov & Saalfrank 2006
            case 'fermi'
                util.disp ('Relaxation rates from Fermi''s golden rule')
                util.disp (['relaxation rate = ' num2str(control.relax.rate)])
                util.disp (['==> upper state = ' int2str(control.relax.upper)])
                util.disp (['==> lower state = ' int2str(control.relax.lower)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.relax.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.relax.lower) '), '] ;

                ratios = abs( triu(tise.sbc,1) ).^2 / abs( tise.sbc(control.relax.lower+1,control.relax.upper+1) )^2;
                
                if kBT == 0 % for zero temperature: only downward relaxation
                    for k=1:dim % upper right triangle of Gamma matrix
                        for l=k+1:dim % l \geq k
                            Gam(k,l)=ratios(k,l) ...
                                *w(control.relax.upper+1,control.relax.lower+1) / w(l,k) ...
                                * control.relax.rate;
                        end
                    end
                    
                else % for finite temperature
                    for k=1:dim % upper right triangle: downward relaxation
                        for l=k+1:dim % l \geq k
                            Gam(k,l) = ratios(k,l) ...
                                * w(control.relax.upper+1,control.relax.lower+1) / w(l,k)...
                                * (1-exp(-w(control.relax.upper+1,control.relax.lower+1)/kBT)) ...
                                / (1-exp(-w(l,k)/kBT)) ...
                                * control.relax.rate;
                            % upward transitions from detailed balance, see Eq. (4)
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from Einstein's spontaneous emission
            case 'einstein'
                util.disp ('Relaxation rates from Einstein''s spontaneous emission')
                util.disp (['relaxation rate = ' num2str(control.relax.rate)])
                util.disp (['==> upper state = ' int2str(control.relax.upper)])
                util.disp (['==> lower state = ' int2str(control.relax.lower)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.relax.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.relax.lower) '), '] ;
                
                ratios = abs( triu(tise.d_x,1) ).^2 / abs( tise.d_x(control.relax.lower+1,control.relax.upper+1) )^2;
                
                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l) = ratios(k,l) ...
                            * (w(l,k) / w(control.relax.upper+1,control.relax.lower+1))^3 ...
                            * control.relax.rate;
                        % upward transitions from detailed balance, see Eq. (4)
                        Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                    end
                end
                
            % Constant relaxation rates (for all downward transitions)
            case 'constant'
                util.disp ('Constant relaxation rates')
                util.disp (['relaxation rate = ' num2str(control.relax.rate)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (const.), '] ;

                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l)=control.relax.rate;
                        if kBT>0 % upward transitions from detailed balance
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from data file
            case 'datafile'
                util.disp ('Reading relaxation rates from file')
                bilinear.title = [bilinear.title '\Gamma from file, '] ;

                data = load ('relax');
                Gam = data.relax;
                
            case 'none'
                util.disp ('Not accounting for relaxation rates')

            otherwise
                util.error ('Invalid choice of relaxation model')
                
        end
        
        if ~strcmpi(control.relax.model,'none')
            
            % Construct total dephasing rate, i. e. gamma matrix from Eq. (5)
            GamSummed = sum(Gam, 1);
            gamma.r = zeros(dim);
            for k=1:dim
                for l=1:dim
                    gamma.r(k,l)=0.5*(GamSummed(k)+GamSummed(l));
                end
            end
            
            % population gain, similar to Eq. (A2), but for columnwise order
            for l=0:dim-1
                ll = 1 + l*dim + l;
                for k=0:dim-1
                    kk = 1 + k*dim + k;
                    bilinear.A(ll, kk) = bilinear.A(ll, kk) + Gam(l+1,k+1);
                end
            end
            
            % Include population loss and total dephasing in matrix A
            bilinear.A = bilinear.A + diag( reshape(-gamma.r, dim^2, 1) );
            
            % Find extrema of gamma.r
            min_gr = abs(gamma.r(1,2)); min_lo=1-1; min_up=2-1;
            max_gr = abs(gamma.r(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim % upper right triangle with diag: downward relaxation
                for l=k:dim % l \gt k
                    abs_gr = abs(gamma.r(k,l));
                    if abs_gr < min_gr
                        min_gr = abs_gr;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gr > max_gr
                        max_gr = abs_gr;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end

            util.disp (['min. relax. dephasing = ' num2str(min_gr)])
            util.disp (['==> upper state       = ' int2str(min_up)])
            util.disp (['==> lower state       = ' int2str(min_lo)])
            util.disp (['==> Bohr frequency    = ' num2str(w(min_up+1,min_lo+1))])
            util.disp (['max. relax. dephasing = ' num2str(max_gr)])
            util.disp (['==> upper state       = ' int2str(max_up)])
            util.disp (['==> lower state       = ' int2str(max_lo)])
            util.disp (['==> Bohr frequency    = ' num2str(w(max_up+1,max_lo+1))])
                 
        end
        util.disp (' ')      
        
        % Lindblad operator for pure dephasing
        if ~isfield(control,'depha')
            control.depha = [];
        end
        if ~isfield(control.depha,'model')
            control.depha.model = 'none';
        end
        
        switch lower (control.depha.model)
            
            % Dephasing rates from stochastic Gaussian model (quadratic energy gap dependence)
            case 'gauss'
                util.disp ('Dephasing rates from stochastic Gaussian model')
                util.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                util.disp (['==> upper state     = ' int2str(control.depha.upper)])
                util.disp (['==> lower state     = ' int2str(control.depha.lower)])
                util.disp (['==> Bohr frequency  = ' num2str(w(control.depha.upper+1,control.depha.lower+1))])
                bilinear.title = [bilinear.title '\Gamma^* =' num2str(control.depha.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.depha.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.depha.lower) '), '] ;

                gamma.d = (w/w(control.depha.upper+1,control.depha.lower+1)).^2 * control.depha.rate;
                
                % Constant dephasing rates
            case 'constant'
                util.disp ('Constant pure dephasing rates')
                util.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                bilinear.title = [bilinear.title '\Gamma^* =' num2str(control.depha.rate)];
                bilinear.title = [bilinear.title ' (const.), '] ;

                gamma.d = ones(size(w)) * control.depha.rate;
                
                % Read dephasing rates from data file
            case 'datafile'
                util.disp ('Reading pure dephasing data from file')
                bilinear.title = [bilinear.title '\Gamma^* from file, '] ;

                data = load ('depha');
                gamma.d = data.depha;
                
            case 'none'
                util.disp ('Not accounting for pure dephasing rates')
                
            otherwise
                util.error ('Invalid choice of pure dephasing model')
                
        end
        
        if ~strcmpi(control.depha.model,'none')
            
            % Contribution to A-matrix
            bilinear.A = bilinear.A + diag( reshape(-gamma.d, dim^2, 1) );
            
            % Find extrema
            min_gd = abs(gamma.d(1,2)); min_lo=1-1; min_up=2-1;
            max_gd = abs(gamma.d(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim-1 % upper right triangle: downward relaxation
                for l=k+1:dim % l \geq k
                    abs_gd = abs(gamma.d(k,l));
                    if abs_gd < min_gd
                        min_gd = abs_gd;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gd > max_gd
                        max_gd = abs_gd;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end
            
            util.disp (['min. pure dephasing = ' num2str(min_gd)])
            util.disp (['==> upper state     = ' int2str(min_up)])
            util.disp (['==> lower state     = ' int2str(min_lo)])
            util.disp (['==> Bohr frequency  = ' num2str(w(min_up+1,min_lo+1))])
            util.disp (['max. pure dephasing = ' num2str(max_gd)])
            util.disp (['==> upper state     = ' int2str(max_up)])
            util.disp (['==> lower state     = ' int2str(max_lo)])
            util.disp (['==> Bohr frequency  = ' num2str(w(max_up+1,max_lo+1))])
            
        end
        util.disp (' ')
        
        % Set up matrix N, similar to  Eq. (A3), but for columnwise order
        if isfield(tise,'d_x')
            bilinear.N{1}=zeros(dim^2);
        end
        if isfield(tise,'d_y')
            bilinear.N{2}=zeros(dim^2);
        end
        
        for l=0:dim-1
            for m=0:dim-1
                index_lm = 1 + l + m*dim;
                
                for k=0:dim-1
                    index_km = 1 + k + m*dim;
                    index_lk = 1 + l + k*dim;
                    
                    if isfield(tise,'d_x')
                        bilinear.N{1}(index_lm,index_km) = bilinear.N{1}(index_lm,index_km) + tise.d_x(l+1,k+1);
                        bilinear.N{1}(index_lm,index_lk) = bilinear.N{1}(index_lm,index_lk) - tise.d_x(k+1,m+1);
                    end
                    
                    if isfield(tise,'d_y')
                        bilinear.N{2}(index_lm,index_km) = bilinear.N{2}(index_lm,index_km) + tise.d_y(l+1,k+1);
                        bilinear.N{2}(index_lm,index_lk) = bilinear.N{2}(index_lm,index_lk) - tise.d_y(k+1,m+1);
                    end
                end
            end
        end
        
        % Set up B, similar to Eq. (A5), but for columnwise order
        for len=1:length(bilinear.N)
            bilinear.B{len} = bilinear.N{len}*bilinear.x.equilib;
        end
        
        % Choice of observables is given in control.observe.targets
        % Set up C vectors: columnwise ordering of matrices
        for len=1:length(control.observe.targets)
            bilinear.label{len} = tise.lab{control.observe.targets(len)};
            switch tise.obs
                case 'amo'
                    util.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' bilinear.label{len}])
                case 'prj'
                    util.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' bilinear.label{len}])
                otherwise
                    util.error ('Invalid choice of observable for LvNE')
            end
            % Transpose, vectorize, transpose such that C*x gives Tr(O*rho)
            op_mat = tise.mat{control.observe.targets(len)};
            op_mat = op_mat.';
            op_vec = op_mat(:);
            bilinear.C {len} = op_vec';
            bilinear.Q {len} = false; % use Re<c|x> as observable

        end
        % bilinear.C = bilinear.C'; % make C a row cell vector
        util.disp (' ')
        
        % if desired: transform full matrix problem
        % columnwise -> diagonal first (cw2df)
        if isfield (control.lvne,'order') && strcmpi(control.lvne.order,'df')
            U=util.cw2df(dim);
            
            bilinear.A = U*bilinear.A*U';
            
            for len=1:length(bilinear.C)
                bilinear.C{len} = bilinear.C{len}*U';
            end
            
            for len=1:length(bilinear.N)
                bilinear.N{len} = U*bilinear.N{len}*U';
                bilinear.B{len} = U*bilinear.B{len};
            end
            
            bilinear.x.initial = U*bilinear.x.initial;
            bilinear.x.equilib = U*bilinear.x.equilib;
        end
        
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch tise.obs
                case {'amo','prj'}
                    bilinear.y.equilib(len) = real(bilinear.C{len}*bilinear.x.equilib);
                    bilinear.y.initial(len) = real(bilinear.C{len}*bilinear.x.initial);
            end

        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        bilinear.x.initial = bilinear.x.initial - bilinear.x.equilib;
        
    case 'tdse'
        
        util.disp (' ')
        util.disp (' coming from the time-dependent Schoedinger equation   ')
        util.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
        util.disp (' for quantum systems interacting with electrical fields')
        util.disp ('                                                       ')
        util.disp ('    d                                                  ')
        util.disp (' i -- |psi(t)> = ( H  - F(t) mu ) |psi(t)>             ')
        util.disp ('   dt               0                                  ')
        util.disp ('                                                       ')
        util.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
        util.disp (' mu is the Hermitian matrix with dipole moments        ')
        util.disp (' and F(t) is the electric field.                       ')
        util.disp ('                                                       ')
        util.disp ('  <D>(t)= <psi(t)| D |psi(t)>                          ')
        util.disp ('                                                       ')
        util.disp ('                                                       ')
        util.disp ('-------------------------------------------------------')
        util.disp (' ')
        
        % Initial and equilibrium state vectors/densities for propagation
        oct.tdse (real(tise.ham))
        
        % making title string for plots
        bilinear.title = 'TDSE: ' ;
        
        % Set up A; shift the eigenenergies by the ground state energy.
        bilinear.A = -1i*(diag(tise.ham)-tise.ham(1)*eye(dim));
        
        if isfield(tise,'d_x')
            bilinear.N{1}=tise.d_x;
            bilinear.B{1}=bilinear.N{1}*bilinear.x.equilib;
        end
        if isfield(tise,'d_y')
            bilinear.N{2}=tise.d_y;
            bilinear.B{2}=bilinear.N{2}*bilinear.x.equilib;
        end
        
        
        % Set up C vectors or D matrices
        % Choice of observables is given in control.observe.targets
        for len=1:length(control.observe.targets)
            bilinear.label{len} = tise.lab{control.observe.targets(len)};
            switch tise.obs
                case 'amo'
                    util.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' bilinear.label{len}])
                    bilinear.D{len} = tise.mat{control.observe.targets(len)};
                case 'prj'
                    util.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' bilinear.label{len}])
                    bilinear.D{len} = tise.mat{control.observe.targets(len)};
                case 'ovl'
                    util.disp (['Observable ' int2str(len) ': Populations from overlaps with eigenstates: ' bilinear.label{len}])
                    bilinear.C{len} = tise.vec{control.observe.targets(len)}';
                    bilinear.Q{len} = true; % use |<c|x>|^2 as observable
                otherwise
                    util.error ('Invalid choice of observable for TDSE')
            end
        end
        util.disp (' ')
                
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch tise.obs
                case {'amo', 'prj'}
                    bilinear.y.equilib(len) = dot ( bilinear.x.equilib, bilinear.D{len} * bilinear.x.equilib );
                    bilinear.y.initial(len) = dot ( bilinear.x.initial, bilinear.D{len} * bilinear.x.initial );
                case 'ovl'
                    bilinear.y.equilib(len) = abs ( bilinear.C{len}*bilinear.x.equilib )^2;
                    bilinear.y.initial(len) = abs ( bilinear.C{len}*bilinear.x.initial )^2;
            end
        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        bilinear.x.initial = bilinear.x.initial - bilinear.x.equilib;
        
    otherwise
        util.error (['Wrong choice of equation-of-motion type : ' eom])
        
end

% Sparsity of (complex) matrix A, if density below threshold
density = nnz(bilinear.A)/numel(bilinear.A);
util.disp (['Density of matrix A: ' num2str(density)])
if density<0.1
    bilinear.A = sparse(bilinear.A);
    util.disp ('  ==> using sparse storage scheme')
else
    util.disp ('  ==> using full storage scheme')
end

% Sparsity of (real) matrices N, if density below threshold
for d=1:length(bilinear.N)
    density = nnz(bilinear.N{d})/numel(bilinear.N{d});
    util.disp (['Density of matrix N_'  int2str(d) ': ' num2str(density)])
    if density<0.1
        bilinear.N{d} = sparse(bilinear.N{d});
        util.disp ('  ==> using sparse storage scheme')
    else
        util.disp ('  ==> using full storage scheme')
    end
end

% Plot spectrum of A 
oct.spec_A (mfilename)

% Save ABNCD etc matrices to data file (index 0 stands for unbalanced, untruncated)
util.disp(['Saving matrices A, B, N, C|D and densities to file : ' eom '.mat'])
save (eom, 'bilinear');

% Output clock/date/time
util.clock;

end

