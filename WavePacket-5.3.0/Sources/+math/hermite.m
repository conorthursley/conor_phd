%------------------------------------------------------------------------------
%
% Returns the value of the Hermite polynomial H_m(x). MatLab
% does not offer this polynomial except for additional packages,
% and the implementation is relatively straight-forward using the
% three-term recursion
% H_n(x) = 2x H_{n-1}(x) - 2(n-1) H_{n-2}(x)
% and the starting point H_0(x) = 1, H_1(x) = 2x
% This is not the most efficient among the possible approaches, but
% very simple, and the time for calculating the polynomials is
% negligible except for special cases and dumb code.
%
% Note that the returned Hermite polynomials are normalized with respect
% to the exp(-x^2) inner product, i.e.
% int exp(-x^2) H_m(x) H_n(x) dx  = delta_mn
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function retval = hermite(x, m)

if m == 0
    retval = ones(size(x));
elseif m == 1
    retval = 2*x;
else
    % do the full recursion
    order  = 2;
    oldest = ones(size(x));
    old    = 2*x;
    cur    = 2 * x .* old - 2 * oldest;

    while order < m
        oldest = old;
        old    = cur;
        order  = order+1;
        cur    = 2*x.*old - 2*(order-1)*oldest;
    end

    retval = cur;
end

% normalize
retval = retval / sqrt(sqrt(pi) * 2^m * factorial(m));
