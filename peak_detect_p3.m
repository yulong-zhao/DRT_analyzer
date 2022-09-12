function [peaks, bounds]  = peak_detect_p3(x,y)

% New peak detection module for gradient based detection. 

[~,g_sg] = sgolay(3,15); % For 320 data points


p=2; dt = 1;
dy_sg2 = conv(y, factorial(p)/(-dt)^p * g_sg(:,p+1), 'same');



threshold = max(y)*0.04; % 4% threshold.

[PKS,LOCS,W,P] = findpeaks(y,x,'MinPeakDistance',0.2, 'MinPeakHeight',threshold);
[sgPKS,sgLOCS,sgW,sgP] = findpeaks(-dy_sg2,x,'MinPeakDistance',0.1,'MinPeakProminence',max(y)*5.0000e-04);
[inv_sgPKS, inv_sgLOCS, inv_sgW, inv_sgP] = findpeaks(dy_sg2,x,'MinPeakDistance',0.1,'MinPeakProminence',max(y)*5.0000e-04);

sgPKS_y = interp1(x,y,sgLOCS);


sgLOCS(sgPKS_y<threshold) = [];
sgPKS_y(sgPKS_y<threshold) = [];

sgPKS_y(sgLOCS>(LOCS(end)+0.3)) = [];
sgLOCS(sgLOCS>(LOCS(end)+0.3)) = [];



% threshold = max(y)*0.05; % 5% threshold.
% [pks,locs,~,~] = findpeaks(y,x);

% locs = locs(pks>=threshold);
% pks  = pks(pks>=threshold);

% [~,inv_locs,~,~] = findpeaks(-y,x);

loc_zero_i = find(y>threshold/10,1,'first') - 1;
loc_last_i = find(y>threshold/10,1,'last')  + 1;

loc_zero  = x(loc_zero_i);
loc_last  = x(loc_last_i);

inv_locs2 = [loc_zero; inv_sgLOCS; loc_last];
inv_pks_real = interp1(x,y,inv_locs2);

sgPKS_y = interp1(x,y,sgLOCS);



peaks.x = sgLOCS;
peaks.y = sgPKS_y;

bounds.x = inv_locs2;
bounds.y = inv_pks_real;

end