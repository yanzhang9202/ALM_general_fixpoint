%% Add all functions in this directory and subdir to searching path
addpath(genpath('/Users/Yan_Zhang/Documents/MATLAB/Optimization/ALM_general_fixpoint'), '-end');

switch prob_type
    case 'waterfilling'
        makeplot_wf;
        
    case 'mpc'
        makeplot_mpc;
        
end