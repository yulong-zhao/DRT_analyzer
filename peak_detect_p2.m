function [peaks, bounds]  = peak_detect_p2(x,y)

% This is the auto part.

threshold = max(y)*0.05; % 5% threshold.
[pks,locs,~,~] = findpeaks(y,x);

locs = locs(pks>=threshold);
pks  = pks(pks>=threshold);

[~,inv_locs,~,~] = findpeaks(-y,x);

loc_zero_i = find(y>threshold/10,1,'first') - 1;
loc_last_i = find(y>threshold/10,1,'last')  + 1;

loc_zero  = x(loc_zero_i);
loc_last  = x(loc_last_i);

inv_locs2 = [loc_zero; inv_locs; loc_last];
inv_pks_real = interp1(x,y,inv_locs2);



peaks.x = locs;
peaks.y = pks;

bounds.x = inv_locs2;
bounds.y = inv_pks_real;

end