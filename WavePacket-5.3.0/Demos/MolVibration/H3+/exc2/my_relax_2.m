% Copyright (C) 2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.


%----------------------------------------------------------------------
%
% Imaginary time Chebychev propagator for a quantum state vector
% ==============================================================
%
% By introducing an imaginaary time variable, tau = i*t, the
% time-dependent Schrödinger equation can be cast into a the 
% form of a diffusion equation. Then the propagation causes 
% an initially guessed wave function to relax to the lowest 
% eigenfunction, i.e., the quantum mechanical ground state.
%
% R. Kosloff and H. Tal-Ezer, Chem. Phys. Lett. 127(3), 223, (1986)
%
%                   (  de          )        (            )
% exp(-H*tau) = exp (-(--+e   )*tau)  * exp (-alpha*H    )
%                   (  2   min     )        (        norm)
%
% where
%         2    (       de       )
% H     = -- * (H - I*(--+e   ) )
%  norm = de   (       2   min  ) 
%
%          de * tau
% alpha  = --------        
%             2
%
% TAU     imaginary time step
% EMIN    minimum of the Hamiltonian
% DE      range of the Hamiltonian
% ALPHA   dimensionless parameter 
% Hnorm   Hamiltonian normalized to [-1,1]
% H       original Hamiltonian 
% I       unity operator
%
% Then the exponential of the time evolution operator 
% is expanded in a series of complex Chebychev polynomials.
%
%     (             )     N
% exp ( -alpha*H    )  = Sum c (alpha) * psi (-H    )
%     (         norm)    n=0  n             n   norm 
%  
% 
% where the coefficients c_n are modified Bessel functions 
% and where the psi.n are the real Chebychev polynomials
% which are calculated using the recursion 
%
%    psi  = I   (unity)
%       0
%
%    psi  = -H
%       1     norm
%
%    psi  = -2*H     psi    - psi    , n>1
%       n       norm    n-1      n-2
%
%------------------------------------------------------------------------

function my_relax ( step )
global psi hamilt space time gs

persistent my_gs1 my_gs2

% Like Chebychev, but we always apply the operator
% (1-P) H (1-P) instead of H, where P is the projection on the ground state.
% It can be calculated by expanding a wavefunction
% |Psi> = c_0 |GS> + |Rest>,      (1-P) |Psi> = |Rest>
% with c_0 = <GS|Psi>
% => |Rest> = |Psi> - c_0 |GS>




% First step only: Get expansion coefficients; no propagation yet
if step==1
    util.disp (' ')
    util.disp ('************************************************************')
    util.disp ('Numerical propagator: Chebychev polynomials (imaginary time)')
    util.disp ('************************************************************')

    if time.efield.n_pulse > 0
        util.error ('This propagator is only for time independent Hamiltonians (no electric fields)')
    end

    util.disp (['Kosloff number dE*dt / (2*hbar) ' num2str(time.main.alpha)])
    if time.main.alpha <= 10
        util.error ('Kosloff number should not be below 10')
    end

    % Automatic or user-specified truncation of polynomial series
    switch time.propa.order
        case (0)
           
            % Search backward
            for ii = round(10*sqrt(time.main.alpha)) : -1 : round(sqrt(time.main.alpha))
                if abs(besseli(ii,time.main.alpha))>time.propa.precision
                    break
                end
            end
            time.propa.order = ii;
            util.disp (['Automatic truncation if coefficients drop below : ' num2str(time.propa.precision)])
            util.disp (['Number of polynomials required : ' int2str(time.propa.order)])
        otherwise
            util.disp (['Truncating polynomial series after : ' int2str(time.propa.order)])
    end
    util.disp ( ' ' )
     
    % Use MODIFIED Bessel functiosn to get expansion coefficients
    time.main.coeffs = besseli([0:time.propa.order],time.main.alpha);
    time.main.coeffs(2:time.propa.order+1) = 2 * time.main.coeffs(2:time.propa.order+1);

    % Phase factors to compensate for normalization of Hamiltonian
    time.main.expon = exp( -(hamilt.range.delta/2+hamilt.range.min) * time.main.delta );
  
    
    
my_gs1 = gs{1};
my_gs2 = gs{2};

coeff = sum( conj(my_gs1(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs1;

coeff = sum( conj(my_gs2(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs2;

norm2 = sum(abs(psi.dvr.grid_ND{1}(:)).^2 .* space.dvr.weight_ND(:));
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} / sqrt(norm2);
  
else

    %-----------------------------------------------------------
    %  Zero-th Chebyshev polynomial : phi_0 = 1
    %-----------------------------------------------------------
coeff = sum( conj(my_gs1(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs1;
coeff = sum( conj(my_gs2(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs2;
    for m = 1:hamilt.coupling.n_eqs
        cheby0{m} = psi.dvr.grid_ND{m};
        psi.dvr.sum{m} = time.main.coeffs(1) * cheby0{m};
    end
    

    %-----------------------------------------------------------
    %  First Chebychev polynomial phi_1 = - Hnorm
    %-----------------------------------------------------------

    % Projection before and after
coeff = sum( conj(my_gs1(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs1;
coeff = sum( conj(my_gs2(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs2;
    ket.hamilt(0,0,1);
coeff = sum( conj(my_gs1(:)) .* psi.dvr.new_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.new_ND{1} = psi.dvr.new_ND{1} - coeff * my_gs1;
coeff = sum( conj(my_gs2(:)) .* psi.dvr.new_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.new_ND{1} = psi.dvr.new_ND{1} - coeff * my_gs2;
    
    for m = 1:hamilt.coupling.n_eqs
        cheby1{m} = - psi.dvr.new_ND{m};
        psi.dvr.sum{m} = psi.dvr.sum{m} + time.main.coeffs(2) * cheby1{m};
    end
    

    %-----------------------------------------------
    %  Higher Chebychev polynomials (n>1) by recursion:
    %  phi_n = -2*Hnorm phi_{n-1} - phi_{n-2}
    %-----------------------------------------------
    for k=2:time.propa.order

        for m = 1:hamilt.coupling.n_eqs
            psi.dvr.grid_ND{m} = cheby1{m};
        end
coeff = sum( conj(my_gs1(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs1;
coeff = sum( conj(my_gs2(:)) .* psi.dvr.grid_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.grid_ND{1} = psi.dvr.grid_ND{1} - coeff * my_gs2;
        ket.hamilt(0,0,1);
coeff = sum( conj(my_gs1(:)) .* psi.dvr.new_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.new_ND{1} = psi.dvr.new_ND{1} - coeff * my_gs1;
coeff = sum( conj(my_gs2(:)) .* psi.dvr.new_ND{1}(:) .* space.dvr.weight_ND(:) );
psi.dvr.new_ND{1} = psi.dvr.new_ND{1} - coeff * my_gs2;
        for m = 1:hamilt.coupling.n_eqs
            cheby2{m} = - 2 * psi.dvr.new_ND{m} - cheby0{m};
            psi.dvr.sum{m} = psi.dvr.sum{m} + time.main.coeffs(k+1) * cheby2{m};

            cheby0{m} = cheby1{m};
            cheby1{m} = cheby2{m};
        end

    end
        
    %  Multiply wave function with complex phase factor
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = psi.dvr.sum{m} * time.main.expon;
    end
    
    %  Re-normalize wave function
    fac = 0;
    for m = 1:hamilt.coupling.n_eqs
        fac = fac + sum ( abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:) );
    end

    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m}(:) = psi.dvr.grid_ND{m}(:) / sqrt(fac);
    end
    
    % Fake autocorrelation
    time.acf.short = zeros(time.sub.n,1);

end

