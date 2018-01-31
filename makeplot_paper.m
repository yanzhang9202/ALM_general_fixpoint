%% Put all precision result together
clear;
close all;
clc;

w = warning('off', 'all');

% rng(189)
% clrs = rand(4,3);
clrs = [.9, .4, .9;
        .7, .2, .1;
        .3, .9, .5;
        .8, .3, .8];
mkr = ['o', 's', '*', 'x'];
epsilon = [1:3];
lw = 0.5;

% pb = 'MPC';
pb = 'WF';

figure(101)  % optimiality
hold on
for i = 1 : 3
   filename = ['data/', pb, '_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
%    filename = ['data/MPC_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
   load(filename);   
   ho(i) = plot(fval_avg-fstar, 'LineStyle', '-', 'LineWidth', lw, 'Color', clrs(i,:));
   lg = length(fval_avg);
   plot(1:round(lg/10):lg, fval_avg(1:round(lg/10):lg)-fstar, 'LineStyle', 'none', ...
       'Color', clrs(i,:), 'Marker', mkr(i))
   plot(bd_upper(end)*ones(1,length(fval_avg)), 'LineStyle', '--',...
       'Color', clrs(i,:))
   plot(bd_lower(end)*ones(1,length(fval_avg)), 'LineStyle', '-.',...
       'Color', clrs(i,:))
   set(gca, 'FontSize', 15)
   set(gca, 'Xscale', 'log')
%    set(gca, 'Yscale', 'log') 
end
switch pb
    case 'MPC'
    legend(ho, 'fl - 24, eps - 0.1', 'fl - 28, eps - 0.01', 'fl - 31, eps - 0.001', 'Location', 'southeast')
    case 'WF'
    legend(ho, 'fl - 18, eps - 0.1', 'fl - 21, eps - 0.01', 'fl - 25, eps - 0.001', 'Location', 'northeast')
end
hold off
% ylim([-0.2, 0.55])
box on;
print([pb, '_opt'], '-depsc');

figure(102)
hold on
for i = 1 : 3
   filename = ['data/', pb, '_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
%    filename = ['data/MPC_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
   load(filename);   
   hf(i) = plot(feas_avg, 'LineStyle', '-', 'LineWidth', lw, 'Color', clrs(i,:));
   lg = length(feas_avg);   
   plot(1:round(lg/10):lg, feas_avg(1:round(lg/10):lg), 'LineStyle', 'none', ...
       'Color', clrs(i,:), 'Marker', mkr(i))   
   plot(bd_feas(end)*ones(1,length(feas_avg)), 'LineStyle', '--',...
       'Color', clrs(i,:))
   set(gca, 'FontSize', 15)
   set(gca, 'Xscale', 'log')
   set(gca, 'Yscale', 'log')
end
switch pb
    case 'MPC'
    legend(ho, 'fl - 24, eps - 0.1', 'fl - 28, eps - 0.01', 'fl - 31, eps - 0.001', 'Location', 'southeast')
    case 'WF'
    legend(ho, 'fl - 18, eps - 0.1', 'fl - 21, eps - 0.01', 'fl - 25, eps - 0.001', 'Location', 'northeast')
end
hold off
ylim([1e-4, 0.5])
box on;
print([pb, '_feas'], '-depsc');

% figure(101)
% legend(ho, 'fl - 24, eps - 0.1', 'fl - 28, eps - 0.01', 'fl - 31, eps - 0.001', 'Location', 'southeast')
% 
% figure(102)
% legend(hf, 'fl - 24, eps - 0.1', 'fl - 28, eps - 0.01', 'fl - 31, eps - 0.001', 'Location', 'southwest')

%%
figure(103)  % optimiality
hold on
for i = 1 : 3
   filename = ['data/', pb, '_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
%    filename = ['data/MPC_ws_eps1e-',num2str(epsilon(i)),'_011718.mat'];
   load(filename);   
   ho(i) = plot(abs(fval_avg-fstar), 'LineStyle', '-', 'LineWidth', lw, 'Color', clrs(i,:));
   lg = length(fval_avg);
   plot(1:round(lg/10):lg, fval_avg(1:round(lg/10):lg)-fstar, 'LineStyle', 'none', ...
       'Color', clrs(i,:), 'Marker', mkr(i))
   plot(ALMparam.epsilon*ones(1,length(fval_avg)), 'LineStyle', '--',...
       'Color', clrs(i,:))
%    plot(bd_lower(end)*ones(1,length(fval_avg)), 'LineStyle', '-.',...
%        'Color', clrs(i,:))
   set(gca, 'FontSize', 15)
   set(gca, 'Xscale', 'log')
   set(gca, 'Yscale', 'log') 
end
switch pb
    case 'MPC'
    legend(ho, 'fl - 24, eps - 0.1', 'fl - 28, eps - 0.01', 'fl - 31, eps - 0.001', 'Location', 'southeast')
    case 'WF'
    legend(ho, 'fl - 18, eps - 0.1', 'fl - 21, eps - 0.01', 'fl - 25, eps - 0.001', 'Location', 'northeast')
end
hold off
% ylim([-0.2, 0.55])
box on;
print([pb, '_opt'], '-depsc');
