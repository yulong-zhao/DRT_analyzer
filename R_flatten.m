function R = R_flatten(R_shaped, cand_info)

aaa = [0,cumsum(cand_info.n_funs_all)];

if(islogical(R_shaped(1)))
    R = false(aaa(end),1);
else
    R = zeros(aaa(end),1);
end




for iii=1:cand_info.n_cand
    R(aaa(iii)+1:aaa(iii+1)) = R_shaped(iii,cand_info.indices{iii});
end

end