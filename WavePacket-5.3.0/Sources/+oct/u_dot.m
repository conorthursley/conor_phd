function u_dot (action,step)
%--------------------------------------------------------------------------
% Derivative of optimal field for propagating BACKWARD/FORWARD
%
% du /dt = - Im ( <z(t) | N A - A N | x(t)+x > ) *s (t)/alpha
%   k                      k       k        e      k         k
%
% see, e.g., Eq. (45) in DOI:10.1063/1.476575
% where one term is missing because u=v was assumed there 
%
% optional variant when optimizing |<c|x>|^2
% du(t)/dt = - Im(<x(t)+x_e|z(t)>*<z(t)|NA-AN|x(t)+x_e>)*s(t)/alpha
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

% Return if second order not rquired
if control.optimal.order<2
    return
end

% So far only one-dimensional control implemented
if control.u.dim>1
    util_error ('multi-dimensional control not yet implemented here!')
end

u = control.u.forward (:,step);
v = control.u.backward(:,step); % in the literatur often called: \bar{u}
x = control.x.forward (:,step) + control.x.equilib;
z = control.x.backward(:,step);

% Calculate matrix elements
braket=zeros(control.u.dim,1);
for d=1:control.u.dim
    matrix = bilinear.N{d}*bilinear.A-bilinear.A*bilinear.N{d};
    matrix = matrix - 1i * (v(d)-u(d)) * bilinear.N{d}^2;
    braket(d) = dot(z,matrix*x);
end

% Optionally multiply with extra factor <x(t)|z(t)>
if isfield(bilinear,'C') && bilinear.Q{control.optimal.terminal}
    switch lower (control.optimal.prefactor)
        case 'initial' % for t=0
            x = control.x.forward (:,1) + control.x.equilib;
            z = control.x.backward(:,1);
        case 'current'
        case 'final' % for t=T
            x = control.x.forward (:,end) + control.x.equilib;
            z = control.x.backward(:,end);
    end
    braket = braket * dot(x,z);
end

% Calculate derivative of control field
derivative = - imag(braket) ...
    .* control.u.shaping(:,step) ...
    ./ control.optimal.alpha;

switch lower(action)
    case ('backward')
        
        % Derivative of control field for propagating BACKWARD
        control.d.backward(:,step) = derivative;
        
        % Mixing with field from previous (forward!) propagation
        if control.optimal.eta~=1
            control.d.backward(:,step) ...
                = control.d.backward(:,step) *    control.optimal.eta ...
                + control.d.forward (:,step) * (1-control.optimal.eta);
        end
        
    case ('forward')
        
        % Derivative of control field for propagating FORWARD
        control.d.forward(:,step) = derivative;
        
        % Mixing with field from previous (backward!) propagation step
        if control.optimal.zeta~=1
            control.d.forward(:,step) ...
                = control.d.forward (:,step) *    control.optimal.zeta ...
                + control.d.backward(:,step) * (1-control.optimal.zeta);
        end
        
    otherwise
        util.error (['Calling u_dot function with wrong "ACTION" keyword: ' action])
end

