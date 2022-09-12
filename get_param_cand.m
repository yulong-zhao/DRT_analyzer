function x_index = get_param_cand(cand, cand_info)

param = (1:cand_info.n_param)';

x_index = cell(cand_info.n_space,1);

R_shaped = R_shaper((1:cand_info.n_space)', cand_info);

bnd2 = 0;
param_count = 0;

if( cand_info.is_same_mu)
    mu_all = reshape(param(1:cand_info.n_peaks),1,[]);
    param_count = cand_info.n_peaks;
end



for i=1:cand_info.n_cand
    
    bnd1 = bnd2 + 1;
    bnd2 = bnd1 - 1 + cand_info.n_funs_all(i);
    
    this_params = param(param_count+1 : param_count+cand_info.n_param_per_cand(i));
    
    
    
    reshaped_param = reshape(this_params, cand(i).n_param, []);
    if(~isempty(reshaped_param))
      %  reshaped_param = [reshaped_param; mu_all(cand_info.indices{i})];
        for ix = bnd1:bnd2
            x_index{ix} = reshaped_param(:,ix-bnd1+1);
            
        end
     %   x_index(bnd1:bnd2)= cand(i).fun(x, reshaped_param);
    end
    
    param_count = param_count+cand(i).n_p_tot;
end



end