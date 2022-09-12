function [Rvals] = recursive_enum(bnb_indices_list, disable_vec, cand_info, cand, Rvals)


if(isempty(bnb_indices_list))
    indices = R_shaper((1:cand_info.n_space)',cand_info);
    R_indices = indices(disable_vec);
    R_indices = R_indices(R_indices~=0);
    
    Rvals = [Rvals, R_indices];
    
else
    
    i = bnb_indices_list(1);
    indices = R_shaper((1:cand_info.n_space)',cand_info);
    j_indices = find(indices(:,i)~=0);
    
    for j = (j_indices(:))'
        
        disable_vec_new = disable_vec;
        
        disable_vec_new(:,i) = true;
        
        disable_vec_new(j,i) = false;
        
        Rvals = recursive_enum(bnb_indices_list(2:end), disable_vec_new, cand_info, cand, Rvals);
    end
end




end