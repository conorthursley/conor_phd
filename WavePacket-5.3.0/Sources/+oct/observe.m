%--------------------------------------------------------------------------
%
% Calculates values of observables
% which can be linear (bilinear.C) and/or quadratic (bilinear.D) 
%
% y = C * x + x * D * x
%
% Depending on variable 'action' the following is achieved:
% 'initial'  : initialization 
% 'forward'  : forward propagation 
% 'backward' : backward propagation 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2012 Burkhard Schmidt, Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

function observe (action, step)

global bilinear control;

switch lower(action)
    
    %% Only upon first call: get/check number of observables
    case 'initial'
        
        % At least one of C or D should exist
        if ~isfield(bilinear,'C') && ~isfield(bilinear,'D')
            util.error ('either cell vector C or D need to exist')
        end
        
        % However, C and D shouldn't both exist (to be generalized later ...)
        if isfield(bilinear,'C') && isfield(bilinear,'D')
            util.error ('cell vectors C and D should not both exist')
        end
        
        % Lengths of cell vectors C or D (the latter should be Hermitian)
        if isfield(bilinear,'C')
            bilinear.len.CD = length(bilinear.C);
        end
        if isfield(bilinear,'D')
            bilinear.len.CD = length(bilinear.D);
            for len=1:bilinear.len.CD
%                 if ~ishermitian (bilinear.D{len})
%                     util.error ('All observable matrices D should be Hermitian')
%                 end
            end
        end
        
        %% Calculate observables from C and/or D from state vector x(t)
    case 'forward'
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        x = control.x.forward (:,step)+control.x.equilib;        
        y = zeros (bilinear.len.CD,1);
        for len = 1:bilinear.len.CD
            
            if isfield(bilinear,'C') % linear observable: <c|x> may be complex
                if ~bilinear.Q{len}
                    y(len) = y(len) + real(bilinear.C{len}*x);
                else
                    y(len) = y(len) +  abs(bilinear.C{len}*x)^2;
                end
            end
            if isfield(bilinear,'D') % quadratic observable: real if D is Hermitian
                y(len) = y(len) + dot( x, bilinear.D{len} * x);
            end
        end
        control.y.forward(:,step) = y;
        
        %% Calculate observables from C and/or D from Lagrange multiplier z(t)
    case 'backward'
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        z = control.x.backward(:,step);
        y = zeros (bilinear.len.CD,1);
        for len = 1:bilinear.len.CD
            if isfield(bilinear,'C') % linear observable: <c|z> may be complex
                if ~bilinear.Q{len}
                    y(len) = y(len) + real(bilinear.C{len}*z);
                else
                    y(len) = y(len) +  abs(bilinear.C{len}*z)^2;
                end
            end
            if isfield(bilinear,'D') % quadratic observable: real if D is Hermitian
                y(len) = y(len) + dot( z, bilinear.D{len} * z);
            end
        end
        control.y.backward(:,step) = y;
        
    otherwise
        util.error('Calling observe.m with invalid keyword')
end

end
