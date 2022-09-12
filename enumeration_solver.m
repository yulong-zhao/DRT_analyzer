sdpv.x = optimvar('x',cand_info.n_param);
sdpv.R = optimvar('R',cand_info.n_space);


sdpv.A = get_kernel_space(x_try, cand, sdpv.x, cand_info);
sdpv.y = sdpv.A*sdpv.R;

%sdpv.int = sum(sum(sdpv.A).*(sdpv.R.^2)'./(0.05+sdpv.R'))/size(sdpv.A,1);

sdpv.err = (sdpv.y - y_try)'*(sdpv.y - y_try);


sdpv.obj = 300*sdpv.err/length(y_try); %0.001*sum(sdpv.R)/length(sdpv.R);


sdpv.const = [0<= sdpv.R; sdpv.R<=1.2; opt.lb <= sdpv.x; sdpv.x <= opt.ub];



[sdpv.cand_for_peak, sdpv.bnb_indices_list] = get_cand_for_peak(cand_info);


sdpv.R_shaped = R_shaper((1:cand_info.n_space)', cand_info);


indices = R_shaper((1:cand_info.n_space)',cand_info);
disable_vec = false(size(indices ));

%%
sdpv.R_vals  = recursive_enum(sdpv.bnb_indices_list, disable_vec, cand_info, cand, []);
%%

sdpv.N_enum = size(sdpv.R_vals,2);

sdpv.results = cell(sdpv.N_enum,1);


for i=1:sdpv.N_enum 

sdpv.warm = warm.all_x0;

sdpv.warm_peak  = sdpv.warm(1:cand_info.n_peaks);

sdpv.x0.x = sdpv.warm;
sdpv.x0.R = interp1(x_try, y_try, sdpv.warm_peak(horzcat(cand_info.indices{:})));


sdpv.x0.R(sdpv.R_vals(:,i)) = 0;

sdpv.R_const = [sdpv.R(sdpv.R_vals(:,i))==0];

sdpv.prob = optimproblem("Objective",sdpv.obj);
sdpv.prob.Constraints.const   = sdpv.const;
sdpv.prob.Constraints.R_const = sdpv.R_const;
sol.opt_val = solve(sdpv.prob, sdpv.x0);

sol.x = sol.opt_val.x;
sol.R = sol.opt_val.R;


sol.A = get_kernel_space(x_try, cand, sol.x , cand_info);
sol.y = sol.A*sol.R;
sol.y_all = sol.A.*sol.R';

sol.err = (sol.y - y_try)'*(sol.y - y_try); 
sol.obj = 200*sol.err/length(y_try) + 0.001*sum(sol.R)/length(sol.R);

sol.enum = sdpv.R_vals(:,i);


%figure; plot(sol.y_all); hold on; plot(y_try)  

sdpv.results{i} = sol;
    
end

[~, sdpv.min_i] = min(cellfun(@(x) x.obj, sdpv.results));
sdpv.best = sdpv.results{sdpv.min_i};


A_opt = sdpv.best.A;
y_hat_opt = sdpv.best.y;
cost_opt  = sdpv.best.obj;
R_opt     = sdpv.best.R;
x_opt     = sdpv.best.x;