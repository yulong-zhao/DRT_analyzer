bnb.indices = R_shaper((1:cand_info.n_space)',cand_info);
bnb.x_index = get_param_cand(cand, cand_info);
bnb.cand_for_peak = cell(cand_info.n_peaks,1);



bnb_disable_vec = false(size(bnb.indices ));


for i_cand = 1:cand_info.n_cand
   for iy = 1:cand_info.n_peaks
       
       if (ismember(iy,cand_info.indices{i_cand}))
       bnb.cand_for_peak{iy} = [bnb.cand_for_peak{iy},i_cand]; 
       end
   end

end

bnb.len_cand_for_peak = cellfun(@length, bnb.cand_for_peak);
assert(all(bnb.len_cand_for_peak>0));
bnb_indices_list = find(cellfun(@length, bnb.cand_for_peak)~=1); % more than one option. 


bnb.opt_fun =  @(param, disable_vec) find_Rbnb( get_kernel_space(x_try, cand, param, cand_info), y_try, 0, cand_info, disable_vec);
bnb.fmincon_fun  = @(x0,disable_vec) call_fmincon(x0, disable_vec, bnb, opt, opti, cand_info); 


bnb.x0_ref = warm.all_x0;
bnb.best.cost = Inf; 
bnb.best.x0 = [];
bnb.best.disable_vec = [];


[cost_last, disable_vec_last, x0_end] = recursive_BNB(bnb.x0_ref, bnb_indices_list, bnb_disable_vec, bnb, x_try, y_try, cand_info, cand);


x_opt = x0_end;




function [cost_last, disable_vec_last, x0_end] = recursive_BNB(x0, bnb_indices_list, disable_vec, bnb, x_try, y_try, cand_info, cand)

for i = (bnb_indices_list(:))'
    x00 = bnb.fmincon_fun(x0, disable_vec); % optimizasyon.
    
    A_opt = get_kernel_space(x_try, cand, x00, cand_info);
    [cost_opt, R_opt, y_hat_opt,res_opt] = find_Rbnb(A_opt, y_try, 0, cand_info, disable_vec);
    
    if(bnb.best.cost < cost_opt)
        cost_last = bnb.best.cost;
        disable_vec_last = bnb.best.disable_vec;
        x0_end = bnb.best.x0;
        break;
    end
    
    
    R_opt_shaped = R_shaper(R_opt, cand_info);
    
    
    
    [~,j_all_ind] = sort(R_opt_shaped(:,i));
    j_indices = find(bnb.indices(:,i)~=0);
    
    j_indices = intersect(j_all_ind(:), j_indices(:));
    
    for j = (j_indices(:))'
        
        disable_vec_new = disable_vec;
        
        disable_vec_new(:,i) = true;
        
        disable_vec_new(j,i) = false;
        
        
        [cost_last, disable_vec_last, x0_end] = recursive_BNB(x0, bnb_indices_list(2:end), disable_vec_new, bnb, x_try, y_try, cand_info, cand);
        
        if(cost_last <  bnb.best.cost) 
            bnb.best.cost = cost_last;
            bnb.best.x0 = x0_end;
            bnb.best.disable_vec = disable_vec_last;
        end
        
    end
    
    
    
end
if(isempty(bnb_indices_list))
    x00 = bnb.fmincon_fun(x0, disable_vec); % optimizasyon.
    
    A_opt = get_kernel_space(x_try, cand, x00, cand_info);
    [cost_opt, R_opt, y_hat_opt,res_opt] = find_Rbnb(A_opt, y_try, 0, cand_info, disable_vec);
    
    R_opt_shaped = R_shaper(R_opt, cand_info);
    
    cost_last = cost_opt;
    disable_vec_last = disable_vec;
    x0_end = x00;
end
end



function x_opt = call_fmincon(x0, disable_vec, bnb, opt, opti, cand_info)

R_indices = bnb.indices(disable_vec);
R_indices = R_indices(R_indices~=0);

x_indices = vertcat(bnb.x_index{R_indices});

active_indices = (1:length(x0))';
active_indices(x_indices) = [];

x0_trim = x0(active_indices); 

lb_trim = opt.lb(active_indices);
ub_trim = opt.ub(active_indices);


opt_fun = @(param) bnb.opt_fun(x_augmenter(param, x0, active_indices), disable_vec);

x_opt_trim = fmincon(opt_fun, x0_trim, opt.A, opt.b, opt.Aeq, opt.beq, lb_trim, ub_trim, opt.nonlcon, opti.fmincon);

x_opt = x_augmenter(x_opt_trim, x0, active_indices);


end


function x_aug = x_augmenter(x0_trim, x0_ref, active_indices)

x_aug = x0_ref;
x_aug(active_indices) = x0_trim;
end