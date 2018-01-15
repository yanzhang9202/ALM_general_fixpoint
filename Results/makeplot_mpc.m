fig = 1;
% Process the result data
fstar = data.sol.fval;

H = double(data.H);
fval_avg = 1/2*diag(var.x_avg'*H*var.x_avg)/P1; % Scaled back to the actual optimality
fval_iter = 1/2*diag(var.x'*H*var.x)/P1;
clear H

Aeq = double(data.Aeq);
beq = double(data.beq);
vec = Aeq*var.x_avg - repmat(beq,1,size(var.x_avg,2));
feas_avg = sqrt(diag(vec'*vec))/P2;
vec = Aeq*var.x - repmat(beq,1,size(var.x,2));
feas_iter = sqrt(diag(vec'*vec))/P2;
clear Aeq beq vec

E = alpha*ALMparam.epsilon;
iter = linspace(1,Ko,Ko);        
bd_upper = (C(1)./iter+E)/P1;   % Scaled back to the normal bound
bd_lower = -(C(2)./iter+E)/P1;
bd_feas = (C(3)./iter+E)/P1;
% Plot for the primal optimality
lw = 1;
figure(fig)
hold on;
h(1)=plot(fval_avg-fstar, 'r-', 'LineWidth', lw);     % Of ergodic average
h(2)=plot(fval_iter-fstar, 'm-', 'LineWidth', lw);    % Of direct iterate
h(3)=plot(bd_upper, 'b-', 'LineWidth', lw);   % Upper bound
h(4)=plot(bd_lower, 'g-', 'LineWidth', lw);   % Lower bound
legend(h, 'Ergodic average', 'Iterate', 'Upper Bound', 'Lower Bound');
hold off;
box on;
% set(gca, 'Yscale', 'log')
set(gca, 'FontSize', 15);
% create a new pair of axes inside current figure
axes('position',[.35 .175 .45 .3])
box on % put box around new pair of axes
indexOfInterest = [(Ko-10) : Ko]; % range of t near perturbation
hold on;
plot(indexOfInterest, fval_iter(indexOfInterest) - fstar, 'm-'); % plot on new axes
plot(indexOfInterest, fval_avg(indexOfInterest) - fstar, 'r-'); % plot on new axes
plot(indexOfInterest, bd_upper(indexOfInterest), 'b-'); % plot on new axes
plot(indexOfInterest, bd_lower(indexOfInterest), 'g-'); % plot on new axes
% set(gca, 'Yscale', 'log')
set(gca, 'FontSize', 15)
axis tight
hold off;
fig = fig + 1;
h = [];
% Plot for the primal feasiblity  
figure(fig)
hold on;
h(1) = plot(feas_avg, 'r-', 'LineWidth', lw);
h(2) = plot(feas_iter, 'm-', 'LineWidth', lw);
h(3) = plot(bd_feas, 'b-', 'LineWidth', lw);
legend(h, 'Ergodic average', 'Iterate', 'Upper Bound');
hold off;
box on;
set(gca, 'Yscale', 'log')
set(gca, 'FontSize', 15)
% create a new pair of axes inside current figure
axes('position',[.35 .35 .45 .3])
box on % put box around new pair of axes
indexOfInterest = [(Ko-10) : Ko]; % range of t near perturbation
hold on;
plot(indexOfInterest, feas_iter(indexOfInterest), 'm-'); % plot on new axes
plot(indexOfInterest, feas_avg(indexOfInterest), 'r-'); % plot on new axes
plot(indexOfInterest, bd_feas(indexOfInterest), 'b-'); % plot on new axes
xlim([indexOfInterest(1),indexOfInterest(end)])
set(gca, 'Yscale', 'log')
set(gca, 'FontSize', 15)
hold off;      