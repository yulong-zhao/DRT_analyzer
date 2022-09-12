function A = get_kernel_space(x, cand, param, cand_info)

is_same_mu = cand_info.is_same_mu;
n_space = cand_info.n_space;


if(isnumeric(param))
    A = zeros(length(x), n_space);
else
    A = [];
end



bnd2 = 0;
param_count = 0;
if(is_same_mu)
    mu_all = reshape(param(1:cand_info.n_peaks),1,[]);
    param_count = cand_info.n_peaks;
end



for i=1:cand_info.n_cand
    
    bnd1 = bnd2 + 1;
    bnd2 = bnd1 - 1 + cand(i).n_funs;
    
    this_params = param(param_count+1 : param_count+cand(i).n_p_tot);
    
    
    
    reshaped_param = reshape(this_params, cand(i).n_param, []);
    if(~isempty(reshaped_param))
        reshaped_param = [reshaped_param; mu_all(cand_info.indices{i})];
        
        if(isnumeric(reshaped_param))
        A(:, bnd1:bnd2)= cand(i).fun(x, reshaped_param);
        else
            for bn = 1:(bnd2-bnd1+1)
                A = [A, cand(i).fun(x, reshaped_param(:,bn))];
            end
        end
    end
    
    param_count = param_count+cand(i).n_p_tot;
end

end