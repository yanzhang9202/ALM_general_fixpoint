%         if ~isempty(find(GPMflag==2))
%             fprintf('Warning: There are GPM steps not successfully solved!\n')
%         end
%         if ALMparam.exitflag == 1
%             fprintf(['Info: ALM succeeds after ', num2str(k), ' iterations!\n'])
%         end
fig = 1;
% Process the result data
a = double(data.a);
fstar = sum(-log(1./(a + data.sol.x)));

fval_avg = sum(-log(1./(repmat(a,1,Ko) + var.x_avg)),1);
fval_iter = sum(-log(1./(repmat(a,1,Ko) + double(var.x))),1);

feas_avg = abs(ones(1,data.N)*var.x_avg - double(data.pw));
% feas_avg = feas_avg/P2;
feas_iter = abs(ones(1,data.N)*double(var.x) - double(data.pw));
% feas_iter = feas_iter/P2;

clear a

E = alpha*ALMparam.epsilon;
iter = linspace(1,Ko,Ko);        
bd_upper = C(1)./iter+E;
bd_lower = -(C(2)./iter+E);
bd_feas = C(3)./iter+E;
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
set(gca, 'FontSize', 15);
% create a new pair of axes inside current figure
axes('position',[.35 .175 .45 .3])
box on % put box around new pair of axes
rg = max(10, round(Ko/10));
indexOfInterest = [(Ko-rg) : Ko]; % range of t near perturbation
hold on;
plot(indexOfInterest, fval_iter(indexOfInterest) - fstar, 'm-'); % plot on new axes
plot(indexOfInterest, fval_avg(indexOfInterest) - fstar, 'r-'); % plot on new axes
plot(indexOfInterest, bd_upper(indexOfInterest), 'b-'); % plot on new axes
plot(indexOfInterest, bd_lower(indexOfInterest), 'g-'); % plot on new axes
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
rg = max(10, round(Ko/10));
indexOfInterest = [(Ko-rg) : Ko]; % range of t near perturbation
hold on;
plot(indexOfInterest, feas_iter(indexOfInterest), 'm-'); % plot on new axes
plot(indexOfInterest, feas_avg(indexOfInterest), 'r-'); % plot on new axes
plot(indexOfInterest, bd_feas(indexOfInterest), 'b-'); % plot on new axes
xlim([indexOfInterest(1),indexOfInterest(end)])
set(gca, 'Yscale', 'log')
set(gca, 'FontSize', 15)
hold off;      
