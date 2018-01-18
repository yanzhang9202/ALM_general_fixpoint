filename = 'MPC_ws_eps1e-3_011718.mat';
load(filename)

H = double(data.H);
fval_avg = zeros(Ko,1);
fval_iter = zeros(Ko,1);
for i = 1 : Ko
   fval_avg(i) = 1/2*var.x_avg(:,i)'*H*var.x_avg(:,i)/P1;
   fval_iter(i) = 1/2*var.x(:,i)'*H*var.x(:,i)/P1;
end
clear H

Aeq = double(data.Aeq);
beq = double(data.beq);
feas_avg = zeros(Ko,1);
feas_iter = zeros(Ko,1);
for i = 1 : Ko
   feas_avg(i) = norm(Aeq * var.x_avg(:,i) - beq)/P2;
   feas_iter(i) = norm(Aeq * var.x(:,i) - beq)/P2;   
end
clear Aeq beq

E = alpha*ALMparam.epsilon;
iter = linspace(1,Ko,Ko);        
bd_upper = (C(1)./iter+E)/P1;   % Scaled back to the normal bound
bd_lower = -(C(2)./iter+E)/P1;
bd_feas = (C(3)./iter+E)/P1;

clear i

save(filename)