function [ flag ] = check_stopcriter_V2(x, lb, ub, g, delta1, delta2)
flag = 1;   % Inialize with flag = 1
for i = 1 : length(x)
    % For active constraints
    if x(i) == lb(i)    % active at lower bound
        if g(i) < -delta1   % Violate the optimality condition
            flag = 0;
            break;
        else if g(i) < delta1   % Cannot differentiate
                if abs(g(i)) >  delta2
                    flag = 0;
                    break
                end
            end
        end
        
    else if x(i) == ub(i)    % active at upper bound
            if g(i) > delta1   % Violate the optimality condition
                flag = 0;
                break;
            else if g(i) > -delta1   % Cannot differentiate
                    if abs(g(i)) >  delta2
                        flag = 0;
                        break
                    end
                end
            end
            
        else if abs(g(i)) > delta2  % For inactive constraints 
                flag = 0;
                break
             end
        end
    end
end
end