data = obj;
H = double(data.H);
f = double(data.f);
Aeq = double(data.Aeq);
beq = double(data.beq);
ub = double(data.ub);
lb = double(data.lb);
% Solve the fixed point problem exactly using Matlab function
options = optimoptions('quadprog',...
'Algorithm','interior-point-convex','Display','off');
[sol.x, sol.fval, sol.flagexit, sol.output, sol.lambda] = quadprog(...
    double(H), double(f), [], [], double(Aeq), double(beq),...
    double(lb), double(ub), [], options);
if sol.flagexit ~= 1
    error('Formulated problem cannot be exactly solved by Matlab!\n')
end