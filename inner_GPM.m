function [sol, err, flag, k] = inner_GPM(lambda, data, GPMparam, FPparam)
global verbose prob_type gpmverbose
stopcriter = GPMparam.stopcriter;  % 1 - stop by criteria, 2 - stop by iter number
iter_max = GPMparam.iter_max;
epsilon = GPMparam.epsilon;
alpha = GPMparam.alpha;
sol_type = GPMparam.sol_type;

delta1 = GPMparam.delta1;
delta2 = GPMparam.delta2;
delta3 = GPMparam.delta3;

% Initialize a starting point
x = fi(zeros(data.N, 1), FPparam.T, FPparam.F);   % Need to be feasible
x_prev = x;
err = 0;
if sol_type == 2
    x_avg = double(x);
end
flag = 1;   % 1 - success, 2 - fail

% if strcmp(prob_type, 'mpc')
%     M3 = data.M3;
%     M4 = data.M4;
%     M5 = data.M5;
% end

switch prob_type
    case 'waterfilling'
        M3 = data.M3;
        M4 = data.M4;
        
    case 'mpc'
        M3 = data.M3;
        M4 = data.M4;
        M5 = data.M5;
end
   

for k = 1 : iter_max
    x_prev = x;
    
    switch prob_type
        case 'waterfilling'
%             g = calc_grad(x, data, lambda);
%             x = x_prev - alpha*g;
            x = x_prev + M3 * (1./(data.a + x_prev)) - M3 * lambda ...
                - M4 * sum(x_prev) + M4;
            
        case 'mpc'
            x = x_prev - M3*x_prev - M4*lambda + M5;
        otherwise
            error('Error:In GPM, undefined step!');
    end
    
    x = calc_proj(x, data.lb, data.ub);
    if sol_type == 2
        x_avg = (x_avg*k + double(x))/(k+1);
    end
    % Check stopping criteria
    if stopcriter == 1
        if norm(x - x_prev) < epsilon
            flag = 1;
            break;
        end
    end
    
    if stopcriter == 3
        switch prob_type
            case 'waterfilling'
                if sol_type == 1
                    g = - M3 * (1./(data.a + x)) + M3 * lambda ...
                + M4 * sum(x) - M4;
%                     ind_lb = find(x==data.lb);
%                     ind_ub = find(x==data.ub);
                else if sol_type == 2
                        g = - M3 * (1./(data.a + fi(x_avg, FPparam.T, FPparam.F)))...
                            + M3 * lambda + M4 * sum(fi(x_avg, FPparam.T, FPparam.F)) - M4;
%                         ind_lb = find(x_avg==data.lb);
%                         ind_ub = find(x_avg==data.ub);
                    end
                end
                
            case 'mpc'
                if sol_type == 1
                    g = M3*x + M4*lambda - M5;
%                     ind_lb = find(x==data.lb);
%                     ind_ub = find(x==data.ub);
                else if sol_type == 2
                        g = M3*fi(x_avg, FPparam.T, FPparam.F) + M4*lambda - M5;
%                         ind_lb = find(x_avg==data.lb);
%                         ind_ub = find(x_avg==data.ub);
                    end
                end
        end
        
%         exitflag = check_stopcriter(double(g), ind_lb, ind_ub, double(delta1), double(delta2));
        
        if sol_type == 1
            exitflag = check_stopcriter_V2(double(x), double(data.lb), ...
                double(data.ub), double(g), double(delta1), double(delta2));
        else if sol_type == 2
                exitflag = check_stopcriter_V2(x_avg, double(data.lb), ...
                double(data.ub), double(g), double(delta1), double(delta2));
            end
        end
        
        if exitflag
            if gpmverbose
                fprintf(['GPM stops after ', num2str(k), ' iterations!\n']);
            end
            break;
        end
    end
    
    if stopcriter == 4
        switch prob_type
            case 'waterfilling'
                if sol_type == 1
                    g = - M3 * (1./(data.a + x)) + M3 * lambda ...
                + M4 * sum(x) - M4;
%                     ind_lb = find(x==data.lb);
%                     ind_ub = find(x==data.ub);
                else if sol_type == 2
                        g = - M3 * (1./(data.a + fi(x_avg, FPparam.T, FPparam.F)))...
                            + M3 * lambda + M4 * sum(fi(x_avg, FPparam.T, FPparam.F)) - M4;
%                         ind_lb = find(x_avg==data.lb);
%                         ind_ub = find(x_avg==data.ub);
                    end
                end
                
            case 'mpc'
                if sol_type == 1
                    g = M3*x + M4*lambda - M5;
%                     ind_lb = find(x==data.lb);
%                     ind_ub = find(x==data.ub);
                else if sol_type == 2
                        g = M3*fi(x_avg, FPparam.T, FPparam.F) + M4*lambda - M5;
%                         ind_lb = find(x_avg==data.lb);
%                         ind_ub = find(x_avg==data.ub);
                    end
                end
        end        
    % Unfinished, check the norm of the gradient smaller than delta3            
    end
    
end

if stopcriter == 1
    err = norm(x - x_prev);
    if err > epsilon
        flag = 2;
    end
end

if sol_type == 1
    sol = x;
else if sol_type == 2
        sol = fi(x_avg, FPparam.T, FPparam.F);
    else
        error('Error: Undefined solution type in GPM!')
    end
end

if verbose
    if flag == 2
        fprintf('Warning:GPM may not return optimal solution!\n');
    end
end

end