%----------------------------------------------------------------------
%
% Imaginary time Chebychev propagator for a quantum state vector
% ==============================================================
%
% By introducing an imaginary time variable, tau = i*t, the
% time-dependent Schrödinger equation can be cast into the 
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
% is expanded in a series of real Chebychev polynomials.
%
%     (             )     N
% exp ( -alpha*H    )  = Sum c (alpha) * phi (-H    )
%     (         norm)    n=0  n             n   norm 
%  
% 
% where the coefficients c_n are modified Bessel functions 
% and where the psi_n are the real Chebychev polynomials
% which are calculated using the recursion 
%
%    phi  = I   (unity)
%       0
%
%    phi  = -H
%       1     norm
%
%    phi  = -2*H     phi    - phi    , n>1
%       n       norm    n-1      n-2
%
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function cheby_imag ( step )
global psi hamilt space time

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
        util.error ('Kosloff number should not be below 10, then it becomes inefficient.')
    end

    % Automatic or user-specified truncation of polynomial series
    if ~isfield (time.propa,'order')
        time.propa.order = 0;
    end

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
     
    % Use MODIFIED Bessel functions to get expansion coefficients
    time.main.coeffs = besseli([0:time.propa.order],time.main.alpha);
    time.main.coeffs(2:time.propa.order+1) = 2 * time.main.coeffs(2:time.propa.order+1);

    % Factors to compensate for normalization of Hamiltonian
    time.main.expon = exp( -(hamilt.range.delta/2+hamilt.range.min) * time.main.delta );

	% Typical values for the wave function are of order 1 (give or take a few
	% orders of magnitude). The typical values of the Chebychev polynomials that
	% we are juggling around range from about 1 to 1/time.main.expon. Since the
	% largest absolute value that can be represented by doubles is around
	% 10^320, we have to check that we do not play around with too large
	% numbers. I use a very generous offset here, because numerical errors seem
	% to appear already earlier.
    if time.main.expon < 1e-200
        util.error('Time step is too large. Chebychev expansion might not converge.');
    end
  
else

    %-----------------------------------------------------------
    %  Zero-th Chebyshev polynomial : phi_0 = 1
    %-----------------------------------------------------------
    for m = 1:hamilt.coupling.n_eqs
        cheby0{m} = psi.dvr.grid_ND{m};
        psi.dvr.sum{m} = time.main.coeffs(1) * cheby0{m};
    end
    

    %-----------------------------------------------------------
    %  First Chebychev polynomial phi_1 = - Hnorm
    %-----------------------------------------------------------
    ket.hamilt(0,0,1);
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
        ket.hamilt(0,0,1);
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

