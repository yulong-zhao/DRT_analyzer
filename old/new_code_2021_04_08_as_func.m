function [] = new_code_2021_04_08_as_func(i_batch, my_save_name)
myData_name = 'data_0_100_info5_cleaned2.mat'; %'data_0_100_info5.mat';

myData = load(myData_name);

is_cleaning = contains(myData_name, 'clean');

data = [myData.t{i_batch}, myData.g_DRT{i_batch}];

%clear find_R quadprog_lasso yalmip_intquad;
% Let's separate fitting linearly part.
%load('data_volkan.mat');
data_wrangling;

i_try = 1;
rng default; % For reproducibility
addpath('cbrewer/cbrewer');

cand_info.mode = "full";


locs_try = [];
%peak_finding;

%[const_loc_temp, const_lb_ub_temp]= peak_detect(x,y);
pause(0.5);
const_loc_temp = myData.peak_sel{i_batch}.locs;
const_lb_ub_temp = myData.peak_sel{i_batch}.lb_ub;

cand_info.n_peaks = length(const_loc_temp);
cand_info.indices = myData.peak_sel{i_batch}.indices;


cand_info.n_funs_all = cellfun(@length, myData.peak_sel{i_batch}.indices); % Number of functions from candidates. 
cand_info.N_PSO_mat = 1500;
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


% n_ga = n_param + n_funs; % Add n_funs number of binary num. 
% 
% warm.alls_x0 = {};
% 
% for i_cand=1:length(cand)
%     warm.alls_x0{i_cand} = cand(i_cand).warm.x0(:);
%     
%     
% end
% warm.all_x0 = vertcat(warm.alls_x0{which_funs});

%%
%create_particles;


%%






%const.lb_ub(2:2:16,:) = const_lb_ub_temp;
% const.lb_ub(2:2:16,:) = const_loc_temp[-3.9, -3.6;
%                          -3.6, -3.14;
%                          -3.14, -2.5;
%                          -2.9,  -2.1;
%                          -2.2,  -1;
%                          -1.1,  -0.6;
%                          -0.6,   2;
%                          -0.6,   2];


%if(is_cleaning)
    warm.all_x0 = myData.x_opt_for_warm{i_batch};
%else
    warm.all_x0 = mean(warm.lb_ub,2);
    warm.all_x0(1:cand_info.n_peaks) = const_loc_temp; %[-3.7941; -3.4075; -2.7919;  -2.5; -1.5752; -0.5205;    -0.0797;    0.3931];   
%end






PSO_InitialSwarmMatrix = warm.all_x0(:)';



opt.A = []; %temp_A;
opt.b = []; %temp_b;
opt.Aeq = [];
opt.beq = [];
opt.lb = const.lb_ub(:,1);
opt.ub = const.lb_ub(:,2);
opt.nonlcon = [];


opti_settings;

create_random_points_for_ms;

%get_kernel_space(x_try,cand,mean(const.lb_ub,2), cand_info)



%%
opt_fun = @(param) find_R( get_kernel_space(x_try, cand, param, cand_info), y_try, 0, cand_info);
opt_fun_res = @(param) find_Res( get_kernel_space(x_try, cand, param, cand_info), y_try, 0, cand_info); % obj_fun(x, y, cand, param);

%x_opt = fminsearch(opt_fun, warm.all_x0, opti.fminsearch);
%x_opt = fmincon(opt_fun, warm.all_x0, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon, opti.fmincon);
%problem = createOptimProblem('fmincon','objective',opt_fun, 'x0', warm.all_x0, 'lb', opt.lb, 'ub', opt.ub, 'options', opti.fmincon); %opt.A, opt.b, opt.Aeq, opt.beq,  opt.nonlcon,
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
ms = MultiStart('Display','iter'); %, 'UseParallel',true
ms.StartPointsToRun = 'bounds';
ms.MaxTime  = 3*60*60; %30 minutes


problem_lsqnonlin = createOptimProblem('lsqnonlin','objective',opt_fun_res, 'x0', warm.all_x0, 'lb', opt.lb, 'ub', opt.ub, 'options', opti.lsqnonlin); %opt.A, opt.b, opt.Aeq, opt.beq,  opt.nonlcon,

tic;
[x_opt, ~] = run(ms, problem_lsqnonlin, InitialPoints.points);
toc








%opt_fun2 = @(param) opt_fun(fmincon(opt_fun, param, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon, opti.fmincon));
%x_opt = surrogateopt(opt_fun, opt.lb, opt.ub);
%x_opt = patternsearch(opt_fun, warm.all_x0, opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, [],  optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf));


% opt.lb = [opt.lb; zeros(n_funs,1)];
% opt.ub = [opt.ub; n_cand*ones(n_funs,1)];
% 
% opt.intcon = (n_param+1):n_ga;

%x_opt = ga(opt_fun,length(warm.all), opt.A, opt.b, opt.Aeq, opt.beq, opt.lb, opt.ub, opt.nonlcon, opt.intcon, opti.ga);
%x_opt = particleswarm(opt_fun, cand_info.n_param, opt.lb, opt.ub, opti.particleswarm);


%%

A_opt = get_kernel_space(x_try, cand, x_opt, cand_info);
[cost_opt, R_opt, y_hat_opt, res_opt] = find_R(A_opt, y_try, 0, cand_info);
fprintf('Value at warm start: %4.4f vs optimal %4.4f\n', opt_fun(warm.all_x0), cost_opt);


%%

% A_opt = get_kernel_space(x_try, cand, x_opt_lsqnonlin2, cand_info);
% [cost_opt, R_opt, y_hat_opt, res_opt] = find_R(A_opt, y_try, 0, cand_info);
% fprintf('Value at warm start: %4.4f vs optimal %4.4f\n', opt_fun(warm.all_x0), cost_opt);

%%
%close all;
f1 = figure; hold on;
%[y_vals, y_vals_all] = total_fun(x,cand, x_opt);
grid on;

plot(x, y,'k');
plot(x_try, y_hat_opt,'r-.');

f2 = figure; hold on;
%[y_vals, y_vals_all] = total_fun(x,cand, x_opt);
grid on;

plot(x, y,'k');
plot(x_try, y_hat_opt,'r-.'); 

for i=1:length(R_opt)
   
    if(R_opt(i)>0.01)
        plot(x_try, A_opt(:,i)*R_opt(i));
        
    end
    
end

saveName = [my_save_name, num2str(i_batch)];
savefig(f1,fullfile('results100',[saveName,'_a','.fig']));
savefig(f2,fullfile('results100',[saveName,'_b','.fig']));
close all;
save(fullfile('results100',[saveName,'.mat']));
end


