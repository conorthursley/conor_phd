%----------------------------------------------------------------------
%
% Real time Chebychev propagator for a quantum state vector
% ===========================================================
%
% Epansion of the real time evolution operator exp (-i H t / hbar )
% in imaginary Chebychev polynomials for time-independent Hamiltonian.
%
% R. Kosloff, J. Phys. Chem. 92(8), 2087, (1988)
%
%                    (    de         )        (              )
% exp(-i*H*dt) = exp (-i*(--+e   )*dt)  * exp (-i*alpha*H    )
%                    (    2   min    )        (          norm)
%
% where
%         2    (       de       )
% H     = -- * (H - I*(--+e   ) )
%  norm = de   (       2   min  ) 
%
%          de * dt
% alpha  = -------        
%             2
%
% DT      time step
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
%     (               )     N
% exp ( -i*alpha*H    )  = Sum a (alpha) * phi (-i*H    )
%     (           norm)    n=0  n             n     norm 
%  
% where the coefficients a_n are Bessel functions and
% where the phi_n are the complex Chebychev polynomials
% which are calculated using the recursion 
%
%    phi  = I   (unity)
%       0
%
%    phi  = -i*H
%       1       norm
%
%    phi  = -2*i*H     phi    + phi    , n>1
%       n         norm    n-1      n-2
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function cheby_real ( step )
global psi hamilt space time

% First step only: Get expansion coefficients; no propagation yet
if step==1
    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ('Numerical propagator: Chebychev polynomials (real time)')
    util.disp ('*******************************************************')

    if time.efield.n_pulse > 0
        util.error ('This propagator is only for time independent Hamiltonians (no electric fields)')
    end

    util.disp (['Kosloff number dE*dt / (2*hbar) ' num2str(time.main.alpha)])
    if time.main.alpha <= 40
        util.error ('Kosloff number should not be below 40, then it becomes inefficient')
    end

    % Automatic or user-specified truncation of polynomial series
    if ~isfield (time.propa,'order')
        time.propa.order = 0;
    end
    
    switch time.propa.order
        case (0)
           
            % Search backward
            for ii = round(2.0*time.main.alpha) : -1 : 0.5*round(time.main.alpha)
                if abs(besselj(ii,time.main.alpha))>time.propa.precision
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
 
    
    
    % Use Bessel functions to get expansion coefficients
    time.main.coeffs = besselj([0:time.propa.order],time.main.alpha);
    time.main.coeffs(2:time.propa.order+1) = 2 * time.main.coeffs(2:time.propa.order+1);

    % Get coefficients (columns) also for time sub steps (rows)
    time.sub.alpha =  [1:time.sub.n]' * time.sub.delta * hamilt.range.delta/2;
   
    % Suddenly(?), BESSELJ  no longer accepts mixtures of row and column vectors
    [NU1,Z1] = meshgrid ([0:time.propa.order],time.sub.alpha);
    time.sub.coeffs  = besselj(NU1,Z1);
    % time.sub.coeffs  = besselj([0:time.propa.order],time.sub.alpha);
    time.sub.coeffs (:,2:time.propa.order+1) = 2 * time.sub.coeffs(:,2:time.propa.order+1);

    % Phase factors to compensate for normalization of Hamiltonian
    time.main.phase = exp( -1i * (hamilt.range.delta/2+hamilt.range.min) *                   time.main.delta );
    time.sub.phase  = exp( -1i * (hamilt.range.delta/2+hamilt.range.min) * [1:time.sub.n]' * time.sub.delta );

else

    %-----------------------------------------------------------
    %  Zero-th Chebyshev polynomial : phi_0 = 1
    %-----------------------------------------------------------
    for m = 1:hamilt.coupling.n_eqs
        cheby0{m} = psi.dvr.grid_ND{m};
        psi.dvr.sum{m} = time.main.coeffs(1) * cheby0{m};
    end
    
    % Autocorrelation
    overlap = 0;
    for m = 1:hamilt.coupling.n_eqs
        overlap = overlap + sum ( conj(psi.dvr.init.ND{m}(:)).*cheby0{m}(:) .* space.dvr.weight_ND(:) );
    end
    time.acf.short = time.sub.coeffs(:,1) * overlap;

    %-----------------------------------------------------------
    %  First Chebychev polynomial phi_1 = -i*Hnorm
    %-----------------------------------------------------------
    ket.hamilt(0,0,1);
    for m = 1:hamilt.coupling.n_eqs
        cheby1{m} = - 1i  * psi.dvr.new_ND{m};
        psi.dvr.sum{m} = psi.dvr.sum{m} + time.main.coeffs(2) * cheby1{m};
    end
    
    % Autocorrelation
    overlap = 0;
    for m = 1:hamilt.coupling.n_eqs
        overlap = overlap + sum ( conj(psi.dvr.init.ND{m}(:)).*cheby1{m}(:) .* space.dvr.weight_ND(:) );
    end
    time.acf.short = time.acf.short + time.sub.coeffs(:,2) * overlap;

    %-----------------------------------------------
    %  Higher Chebychev polynomials (n>1) by recursion:
    %  phi_n = -2*i*Hnorm phi_{n-1} + phi_{n-2}
    %-----------------------------------------------
    for k=2:time.propa.order

        for m = 1:hamilt.coupling.n_eqs
            psi.dvr.grid_ND{m} = cheby1{m};
        end
        ket.hamilt(0,0,1);
        for m = 1:hamilt.coupling.n_eqs
            cheby2{m} = - 2 * 1i * psi.dvr.new_ND{m} + cheby0{m};
            psi.dvr.sum{m} = psi.dvr.sum{m} + time.main.coeffs(k+1) * cheby2{m};

            cheby0{m} = cheby1{m};
            cheby1{m} = cheby2{m};
        end

    % Autocorrelation
        overlap = 0;
        for m = 1:hamilt.coupling.n_eqs
            overlap = overlap + sum ( conj(psi.dvr.init.ND{m}(:)).*cheby2{m}(:) .* space.dvr.weight_ND(:) );
        end
        time.acf.short = time.acf.short + time.sub.coeffs(:,k+1) * overlap;

    end

    %  Multiply wave function with complex phase factor
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = psi.dvr.sum{m} * time.main.phase;
    end
    
    %  Multiply autocorrelation with complex phase factor
    time.acf.short = time.acf.short .* time.sub.phase;

end

