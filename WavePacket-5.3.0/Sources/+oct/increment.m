function delta = increment (iter)
%--------------------------------------------------------------------------
% Increment of total control functional
% see DOI:10.1063/1.1650297
%--------------------------------------------------------------------------
global bilinear control

% Equation (17)
delta = control.optimal.alpha * trapz ...
    (control.t.steps, ...
    (2/control.optimal.zeta-1) * ...
    (control.u.backward - control.u.forward).^2 + ...
    (2/control.optimal.eta-1) * ...
    (control.u.backward - control.u.previous).^2 ./ ...
    control.u.shaping(1,:) );

% Equation (A5)
if isfield(bilinear,'D')
    dxT = control.x.previous(:,end) - control.x.forward(:,end);
    delta = delta + dot( dxT, bilinear.D{control.optimal.terminal} * dxT);
end

% Case (C2), i.e. for D=|c><c|
if isfield(bilinear,'C') && bilinear.Q{control.optimal.terminal}
    dxT = control.x.previous(:,end) - control.x.forward(:,end);
    delta = delta + abs ( bilinear.C{control.optimal.terminal} * dxT )^2;
end

% Output
delta = util.real (delta);
util.disp (['After ' int2str(iter) ' iterations:     DJ = ' num2str(delta)])

end


