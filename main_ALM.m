%% ALM for general problems
clear
close all
% clc
global prob_type ALMparam verbose gpmverbose

%% Initialization
prob_type = 'waterfilling';
% prob_type = 'mpc';
init_problem;       % Now initial the problem with fixed point data
init_param_ALM;

%% Display setting
verbose = 0;
gpmverbose = 0;
do_proj_multiplier = 1;

%% Run ALM
tic
count = 0;
inn_iter = zeros(Ko, 1);
for k = 1 : ALMparam.iter_max
   % Solve inner problem
   switch prob_type
       case 'waterfilling'
           [x_iter, err(k), GPMflag(k), inn_iter(k)] = inner_GPM(lambda_iter,data,GPMparam, FPparam);
       case 'mpc'
%            [var_fp.x(:,k), err(k), GPMflag(k)] = inner_GPM(var_fp.lambda(:,k),data,GPMparam, FPparam);
           [x_iter, err(k), GPMflag(k), inn_iter(k)] = inner_GPM(lambda_iter, data, GPMparam, FPparam);
       otherwise
           error('Error: Undefined solution method for the inner problem!')
   end
   
   var.x(:,k) = double(x_iter);
   var.lambda(:,k) = double(lambda_iter);
   
   % Calculate the ergodic average of the primal and dual var.
   if k == 1
       var.x_avg(:,k) = var.x(:,k);
       var.lambda_avg(:,k) = var.lambda(:,k);
   else
       var.x_avg(:,k) = (var.x_avg(:,k-1)*(k-1)+double(x_iter))/k;
       var.lambda_avg(:,k) = (var.lambda_avg(:,k-1)*(k-1)+...
           double(lambda_iter))/k;
   end   
   % Check stopping condition
   if ALMparam.stopcriter == 1
       if norm(data.A*var.x(:,k)-data.b) < ALMparam.epsilon && GPMflag(k) == 1
           ALMparam.exitflag = 1;
           break;
       end
   end
   
   % Update multiplier
   lambda_prev = lambda_iter;
   switch prob_type
       case 'waterfilling'
%             lambda_iter(:,k+1) = lambda_prev + ALMparam.rho/2*(data.A*x_iter-data.b);
            lambda_iter = lambda_prev + data.M1 * x_iter - data.M2;
            if do_proj_multiplier
                lambda_iter = calc_proj(lambda_iter, data.lb_lambda, data.ub_lambda);
            end
            
       case 'mpc'
            lambda_iter = lambda_prev + data.M1 * x_iter - data.M2;
            if do_proj_multiplier
                lambda_iter = calc_proj(lambda_iter, data.lb_lambda, data.ub_lambda);
            end
   end
   
   count = count + 1;
   
   if mod(k,10) == 1
       fprintf([num2str(k), 'th iteration completed!\n'])
   end
end
toc
%% Show results
makeplot;

%% Save results
filename = 'WF_ws_eps1e-3_011723.mat';
% save(filename, 'ALMparam', 'alpha', 'beta', 'gamma', 'C', 'data', ...
%     'E', 'FPparam', 'Ki', 'Ko', 'P1', 'P2', 'var', 'prob_type');
save(filename);  % For satefy

%% End the code
% end_ALM;