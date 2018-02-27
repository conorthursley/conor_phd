%------------------------------------------------------
%
% Calculate cosine elliptic Mathieu functions of order r, 
% as a function of q and z, using matrix-based methods
% Functions with even r are of period \pi
% Functions with  odd r are of period 2*\pi
% Solutions of Mathieu's equation:  y"+(a-2qcos2z)y = 0
%
% Julio C. Gutierrez-Vega (Monterrey, Mexico)
%
%------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) ?? Julio C. Gutierrez-Vega
%               2009 Burkhard Schmidt
%
% see the README file for license details.

function [c,eigenvalor,coef]=cee(r,q,z)

% INPUT
%   r = 0,1,2,3,.. : integer order
%   q = Mathieu parameter:  real positiv or zero or complex
%   z = independiente variable: in general complex
% OUTPUT
%   c = values of Mathieu function
%   eigenvalor = r-th eigenvalue of Mathieu's equation
%   coef = coefficients of Fourier series


if r<0
   util.error('Order must be semi-positive and real')
end

if mod(r,2)==0
   N=r+50;  %Tamaño de la matriz para calcular los eigenvalores
   M = diag([0:2:2*N].^2,0) + diag(linspace(q,q,N),1) + diag([2*q,linspace(q,q,N-1)],-1);  %Eigen matriz
   
   [A,ets]=eig(M);
   ets=diag(ets);
   pos = r/2 + 1;   %Para el r'esimo orden calcula la posici´on 
   
   [ets1,index]=sort(real(ets));  %Ordena los eigenvalores usando la parte REAL
   eigenvalor = ets(index(pos));  %eigenvalor para el orden r
   coef = transpose(A(:,index(pos)));   %coeficientes de Fourier para el orden r;  arreglados en un renglon

   coef = coef/sqrt(coef(1)^2 + sum(coef.^2)) / sign(coef(1));  %Normalizacion y ajuste de signo
   
   c=0;
   for ii=1:length(coef)
      c=c + coef(ii)*cos(2*(ii-1)*z);
   end
    
else
   N=r+25;  %Tamaño de la matriz para calcular los eigenvalores
   M = diag([1+q,[3:2:2*N+1].^2],0) + diag(linspace(q,q,N),1) + diag(linspace(q,q,N),-1);  %Eigen matriz
   
   [A,ets]=eig(M);
   ets=diag(ets);   
   pos = (r+1)/2;   %Para el r'esimo orden calcula la posici´on 

   [ets1,index]=sort(real(ets));  %Ordena los eigenvalores usando la parte REAL
   eigenvalor = ets(index(pos));  %eigenvalor para el orden r
   coef = transpose(A(:,index(pos)));   %coeficientes de Fourier para el orden r;  arreglados en un renglon

   coef = coef/sqrt(sum(coef.^2)) / sign(coef(1));  %Normalizacion y ajuste de signo
   
   c=0;
   for ii=0:(length(coef)-1)
      c=c + coef(ii+1)*cos((2*ii+1)*z);
   end
end
