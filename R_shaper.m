function R_shaped = R_shaper(R, cand_info)

aaa = [0,cumsum(cand_info.n_funs_all)];

R_shaped = zeros(cand_info.n_cand, cand_info.n_peaks);
for iii=1:cand_info.n_cand
    R_shaped(iii,cand_info.indices{iii}) = R(aaa(iii)+1:aaa(iii+1));
end

end