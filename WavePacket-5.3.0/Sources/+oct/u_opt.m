function u_opt (action,step)
%--------------------------------------------------------------------------
% Optimal field for propagating BACKWARD/FORWARD
%
% u (t) = - Im ( <z(t) | N  | x(t)+x > ) *s (t)/alpha
%  k                      k         e      k         k
%
% optional variant when optimizing |<c|x>|^2
% u(t) = - Im(<x(t)+x_e|z(t)>*<z(t)|N|x(t)+x_e>)*s(t)/alpha
% Evaluating the extra factor for t=T fixes the problem of the
% increments but slows down convergence of the OCT iteration.
%
% Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
%
% optionally mixing with field from previous iteration step
% v(t) =  eta * u(t) + (1- eta) * u~(t) for backward propagation
% u(t) = zeta * u(t) + (1-zeta) * v~(t) for  forward propagation
% see DOI:10.1063/1.1650297
%--------------------------------------------------------------------------
global bilinear control

x = control.x.forward (:,step) + control.x.equilib;
z = control.x.backward(:,step);

% Calculate matrix element
braket = zeros(control.u.dim,1);
for d = 1:control.u.dim
    braket(d) = dot(z,bilinear.N{d}*x);
end

% Optionally multiply with extra factor <x(t)|z(t)>
if isfield(bilinear,'C') && bilinear.Q{control.optimal.terminal}
    switch lower (control.optimal.prefactor)
        case 'initial' % for t=0
            x = control.x.forward (:,1) + control.x.equilib;
            z = control.x.backward(:,1);
        case 'current' % for current t
        case 'final' % for t=T
            x = control.x.forward (:,end) + control.x.equilib;
            z = control.x.backward(:,end);
        otherwise
            util.error(['Wrong choice for prefactor time : ' control.optimal.prefactor])
    end
    braket = braket * dot(x,z);
end

% Calculate optimal control field
field = - imag(braket) ...
            .* control.u.shaping(:,step) ...
            ./ control.optimal.alpha;
        
switch lower(action)
    case ('backward')
        
        % Control field for propagating BACKWARD
        control.u.backward(:,step) = field;
        
        % Mixing with field from (previous!) forward propagation
        if control.optimal.eta~=1
            control.u.backward(:,step) ...
                = control.u.backward(:,step) *    control.optimal.eta ...
                + control.u.forward (:,step) * (1-control.optimal.eta);
        end
        
    case ('forward')
        
        % Control field for propagating FORWARD
        control.u.forward(:,step) = field;
        
        % Mixing with field from (previous!) backward propagation step
        if control.optimal.zeta~=1
            control.u.forward(:,step) ...
                = control.u.forward (:,step) *    control.optimal.zeta ...
                + control.u.backward(:,step) * (1-control.optimal.zeta);
        end
        
    otherwise
        util.error (['Calling u_opt function with wrong "ACTION" keyword: ' action])
end

