function j = functionals (iter)
%--------------------------------------------------------------------------
% Calculate cost, target, and total functionals in optimal control theory
%--------------------------------------------------------------------------
global bilinear control

% Cost functional
j.cost = 0;
for d=1:control.u.dim
    j.cost = j.cost + ...
        control.optimal.alpha(d) ...
        * trapz( ...
        control.t.steps,control.u.forward(d,:).^2 ./ ...
        control.u.shaping(d,:));
end
j.cost = util.real(j.cost);


% Target functional (for linear target C)
if isfield(bilinear,'C')
    j.target = control.y.forward(control.optimal.terminal,end);
    j.target = util.real(j.target);
    if bilinear.Q{control.optimal.terminal}
        j.total = j.target-j.cost;
        util.disp (['After ' int2str(iter) ' iterations (C2): J = ' num2str(j.target) ' - ' num2str(j.cost) ' = ' num2str(j.total)])
    else
        j.total = 2*j.target-j.cost;
        util.disp (['After ' int2str(iter) ' iterations (C1): J = 2*' num2str(j.target) ' - ' num2str(j.cost) ' = ' num2str(j.total)])
    end
end

% Target functional (for quadratic target D)
if isfield(bilinear,'D')
    j.target = control.y.forward(control.optimal.terminal,end);
    j.target = util.real(j.target);
    j.total = j.target-j.cost;
    util.disp (['After ' int2str(iter) ' iterations (D): J = ' num2str(j.target) ' - ' num2str(j.cost) ' = ' num2str(j.total)])
end


