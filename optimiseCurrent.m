function [x_sep, y_sep, y_fit] = optimiseCurrent(appStruct, i_batch)
my_save_name = appStruct.other.SaveName;

data = appStruct.data{i_batch};

data_wrangling;

i_try = 1;
rng default; % For reproducibility
cand_info.mode = "full";

locs_try = [];
%peak_finding;

const_loc_temp = appStruct.peak_info{i_batch}.peaks.x; %myData.peak_sel{i_batch}.locs;
const_lb_ub_temp = appStruct.peak_info{i_batch}.lb_ub;% .x ; % myData.peak_sel{i_batch}.lb_ub;

cand_info.n_peaks = length(const_loc_temp);
cand_info.indices = appStruct.peak_info{i_batch}.indices;


cand_info.n_funs_all = cellfun(@length, appStruct.peak_info{i_batch}.indices); % Number of functions from candidates.
cand_info.is_same_mu = true;

cand_info.n_cand = length(cand_info.n_funs_all);
%cand_info.mu_keep = find(cand_info.n_funs_all>0, 1); % keep mu at this index.

%%

define_candidates;

const.lb_ub = const_lb_ub_temp;
warm.lb_ub = const_lb_ub_temp;

for i_cand=1:cand_info.n_cand
    const.lb_ub = [const.lb_ub; repmat(cand(i_cand).n_limits_eff, cand(i_cand).n_funs,   1)];
    warm.lb_ub = [warm.lb_ub; repmat(cand(i_cand).n_limits_warm_eff, cand(i_cand).n_funs,   1)];
end

%%


warm.all_x0 = mean(warm.lb_ub,2);
warm.all_x0(1:cand_info.n_peaks) = const_loc_temp; %[-3.7941; -3.4075; -2.7919;  -2.5; -1.5752; -0.5205;    -0.0797;    0.3931];

opt.A = []; %temp_A;
opt.b = []; %temp_b;
opt.Aeq = [];
opt.beq = [];
opt.lb = const.lb_ub(:,1);
opt.ub = const.lb_ub(:,2);
opt.nonlcon = [];

create_random_points_for_ms2;

PSO_InitialSwarmMatrix = InitialPoints.all;




opti_settings;



%get_kernel_space(x_try,cand,mean(const.lb_ub,2), cand_info)



%%
opt_fun = @(param) find_R( get_kernel_space(x_try, cand, param(:), cand_info), y_try, 0, cand_info);
opt_fun_res = @(param) find_Res( get_kernel_space(x_try, cand, param, cand_info), y_try, 0, cand_info); % obj_fun(x, y, cand, param);

%x_opt = fminsearch(opt_fun, warm.all_x0, opti.fminsearch);
%x_opt = fmincon(opt_fun, warm.all_x0, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon, opti.fmincon);
%x_opt = lsqnonlin(opt_fun_res, warm.all_x0, opt.lb, opt.ub, opti.lsqnonlin);

%problem = createOptimProblem('lsqnonlin','objective',opt_fun, 'x0', warm.all_x0, 'lb', opt.lb, 'ub', opt.ub); %opt.A, opt.b, opt.Aeq, opt.beq,  opt.nonlcon,

% ms = MultiStart('UseParallel',true,'Display','iter');
% ms.StartPointsToRun = 'bounds';
% ms.MaxTime  = 30*60; %30 minutes
% gs = GlobalSearch; %MultiStart
% gs.NumTrialPoints = 5000;
% gs.StartPointsToRun = 'bounds';
% fprintf('Optimization is started.\n');
%
% [x_opt, ~] = run(gs, problem);

%---
problem_fmincon = createOptimProblem('fmincon','objective',opt_fun, 'x0', warm.all_x0, 'lb', opt.lb, 'ub', opt.ub, 'options', opti.fmincon); %opt.A, opt.b, opt.Aeq, opt.beq,  opt.nonlcon,
problem_lsqnonlin = createOptimProblem('lsqnonlin','objective',opt_fun_res, 'x0', warm.all_x0, 'lb', opt.lb, 'ub', opt.ub, 'options', opti.lsqnonlin); %opt.A, opt.b, opt.Aeq, opt.beq,  opt.nonlcon,

tic;
if strcmpi(appStruct.other.Method, 'globalsearch')
    gs = GlobalSearch('Display','iter');
    gs.StartPointsToRun  = 'bounds';
    gs.MaxTime = appStruct.other.MaxTime; % 1 hour.
    
    [x_opt, ~] = run(gs, problem_fmincon); %,
    
elseif strcmpi(appStruct.other.Method, 'Branch and Bound')
    bnb_optimizer;
    
elseif strcmpi(appStruct.other.Method, 'Enumeration: fmincon')
    enumeration_solver;
    
elseif strcmpi(appStruct.other.Method, 'Particle Swarm Optimisation')
    x_opt = particleswarm(opt_fun,length(warm.all_x0),opt.lb,opt.ub,opti.particleswarm);
    
elseif strcmpi(appStruct.other.Method, 'Genetic Algorithm')
    x_opt = ga(opt_fun,length(warm.all_x0),opt.A, opt.b, opt.Aeq, opt.beq, opt.lb,opt.ub,[], opti.ga);
else
    ms = MultiStart('Display','iter'); %, 'UseParallel',true
    ms.StartPointsToRun = 'bounds';
    ms.UseParallel = 1;
    ms.MaxTime  = appStruct.other.MaxTime; % 3 hours.
    

    [x_opt, ~] = run(ms, problem_lsqnonlin, InitialPoints.points); %,
end
toc
    
A_opt = get_kernel_space(x_try, cand, x_opt, cand_info);
[cost_opt, R_opt, y_hat_opt, res_opt] = find_R(A_opt, y_try, 0, cand_info);

y_fit = y_hat_opt*appStruct.other.norm_factor(i_batch);

for i=1:length(R_opt)
    x_sep = x_try;
    y_sep(:, i) = A_opt(:,i)*R_opt(i)*appStruct.other.norm_factor(i_batch);    
end

%%


%fprintf('Value at warm start: %4.4f vs optimal %4.4f\n', opt_fun(warm.all_x0), cost_opt);


%%

flag_plot_separate = appStruct.other.PlotSeparateFigure;
saveName = [my_save_name, num2str(i_batch)];
save(fullfile(appStruct.abs_path, 'results',[saveName,'.mat']));

if flag_plot_separate

    %close all;
    f1 = figure; hold on;
    %[y_vals, y_vals_all] = total_fun(x,cand, x_opt);
    grid on;
    
    plot(x, y*appStruct.other.norm_factor(i_batch), '-', 'Color', '#0072BD', 'LineWidth', 1);
    plot(x_try, y_fit,'r-.', 'LineWidth', 1);
    xlabel('$\log_{10}(\tau)$', 'fontname', 'times', 'Interpreter','latex')
    ylabel('$G(\tau)$', 'fontname', 'times', 'Interpreter','latex')
    title(['Data set-', num2str(i_batch), '-Fitted'], 'fontname', 'times', 'Interpreter','latex')
    legend('Original data', 'Fitted', 'fontname', 'times')

    f2 = figure; hold on;
    %[y_vals, y_vals_all] = total_fun(x,cand, x_opt);
    grid on;
    
    plot(x, y*appStruct.other.norm_factor(i_batch), '-', 'Color', '#0072BD', 'LineWidth', 1);
    plot(x_try, y_sep,'r-.', 'LineWidth', 1);
    
    xlabel('$\log_{10}(\tau)$', 'fontname', 'times', 'Interpreter','latex')
    ylabel('$G(\tau)$', 'fontname', 'times', 'Interpreter','latex')
    title(['Data set-', num2str(i_batch), '-Separated peaks'], 'fontname', 'times', 'Interpreter','latex')
    legend('Original data', 'Separated peaks', 'fontname', 'times')
    
    
    savefig(f1,fullfile(appStruct.abs_path, 'results',[saveName,'_a','.fig']));
    savefig(f2,fullfile(appStruct.abs_path, 'results',[saveName,'_b','.fig']));
    
%     if(~appStruct.other.ShowPlots)
%         close all;
%     end


end