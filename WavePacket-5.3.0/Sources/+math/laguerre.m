%------------------------------------------------------------------------------
%
% Returns the value of the Laguerre polynomial specified by the
% Rodriguez formula:
%  (a)        x  a  1    ( d  )n  (   -x  n+a  )
% L   (x) =  e  x  ---   (--- )   (  e   x     )
%  n                n!   ( dx )   (            )
%
% Note that these polynomials are not normalized! Also some people use
% e.g.  L_n^{a} * n! or L_n^{a} * (-1)^n.
%
% The calculation can be done by the recurrence relation
% n L_n^{a}(x) = (-x+2n+a-1) L_{n-1}^{a}(x) - (n+a-1) L_{n-2}^{a}(x)
%
% and the starting point
% L_0^{a}(x) = 1,       L_1^{a}(x) = -x + a + 1
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function retval = laguerre(x, n, a)

if a <= -1
    util.error('Generalized Laguerre function is not defined for a <= -1')
end

if any(x < 0)
    util.error('Laguerre functionals usually not used for x < 0; probably an error')
end

if n == 0
    retval = ones(size(x));
elseif n == 1
    retval = -x + a + 1;
else
    oldest = ones(size(x));
    old    = -x + a + 1;
    order  = 2;
    while order <= n
        cur = (-x+2*order+a-1) .* old - (order+a-1).* oldest;
        cur = cur / order;
        oldest = old;
        old = cur;
        order = order+1;
    end
    retval = cur;
end
