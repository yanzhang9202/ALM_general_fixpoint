switch prob_type
    case 'waterfilling'
        epsilon = ALMparam.epsilon;
        bd_lambda = 1;  % By sampling
        main_fpt_wf;
        
    case 'mpc'
        epsilon = ALMparam.epsilon;
%         bd_lambda = 15.5316; % By sampling
        bd_lambda = 15.5316 * P1 / P2;   % The bound of multiplier is improved by scaling
        main_fpt_mpc;
end