function [cand_for_peak, bnb_indices_list] = get_cand_for_peak(cand_info)

cand_for_peak = cell(cand_info.n_peaks,1);


for i_cand = 1:cand_info.n_cand
   for iy = 1:cand_info.n_peaks
       
       if (ismember(iy,cand_info.indices{i_cand}))
       cand_for_peak{iy} = [cand_for_peak{iy},i_cand]; 
       end
   end

end

len_cand_for_peak = cellfun(@length, cand_for_peak);
assert(all(len_cand_for_peak>0));
bnb_indices_list = find(cellfun(@length, cand_for_peak)~=1); % more than one option. 


end