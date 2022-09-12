function peak_info = create_info(x,y)


[possible_peak_locs, peak_lb_ub] = peak_detect(x,y);


n_peaks = length(possible_peak_locs);

myList = {};
for i=1:n_peaks
    myList{i} = ['Peak ', num2str(i)];
end


cand_names = {'Gaussian function', 'trigo function', 'nonsymmetric gaussian'};

indices = {};
tfs = {};

for i=1:length(cand_names)
prString = {'Select the peaks the candidate ', [num2str(i),': ', cand_names{i}], 'is applied: '};


[indx,tf] = listdlg('PromptString',prString, 'ListString',myList);

indices{i} = indx;
tfs{i} = tf;

end

peak_info.locs = possible_peak_locs;
peak_info.lb_ub = peak_lb_ub;
peak_info.indices = indices;
peak_info.tfs = tfs;

end

