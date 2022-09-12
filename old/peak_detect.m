function [possible_peak_locs, peak_lb_ub]  = peak_detect(x,y)
MAX_POINTS = 8;

[pks,locs,w,p] = findpeaks(y,x);

[inv_pks,inv_locs,inv_w,inv_p] = findpeaks(-y,x);

n_pks = length(pks);
n_user_max = max(0, MAX_POINTS-n_pks);


loc_zero_i = find(y>0.001,1,'first');
loc_last_i = find(y>0.001,1,'last');

loc_zero  = x(loc_zero_i);
loc_last  = x(loc_last_i);

inv_locs2 = [loc_zero; inv_locs; loc_last];
inv_pks_real = interp1(x,y,inv_locs2);
ff = figure('Units','Normalized','Position',[0.1 0.1 0.8 0.6]);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

plot(x,y); hold on;
plot(locs,pks,'o');
plot(inv_locs2, inv_pks_real, 'x');








points = [];
user_points = cell(1);


% [x1,y1] = ds2nfu(inv_locs2(1),0);
% [x2,y2] = ds2nfu(inv_locs2(2),1.05);

%flag =

nPatch = length(inv_locs2)-1;
myColors = cbrewer('qual','Set3',nPatch);

for i=1:nPatch
    x1 = inv_locs2(i);
    x2 = inv_locs2(i+1);
    y1 = 0;
    y2 = 1.05;
    
    Xpatch = [x1 x2 x2 x1];
    Ypatch = [y1 y1 y2 y2];
    
    ax = patch(Xpatch,Ypatch,myColors(i,:));%'FaceColor','red','FaceAlpha',.3);
    ax.FaceAlpha = 0.3;
end


title(['Select more break points']);

while true
    temp = ginput(1);
    if(isempty(temp))
        break;
    else
        xline(temp(1), 'LineWidth',2);
        inv_locs2(end+1) = temp(1);
    end
    
end

inv_locs2 = sort(inv_locs2);



% h = rectangle('Position', [x1, y1, x2, y2], ...
%                 'Curvature', 0.2, ...
%                 'FaceColor', [1, 0, 0, 0.7], ...
%                 'EdgeColor', [1, 0, 0, 0.7]);


title(['Select more peaks',newline,'Press enter to finish before maximum of ', num2str(n_user_max), ' points']);

for i=1:n_user_max
    temp = ginput(1);
    
    if(isempty(temp))
        break;
    else
        points = [points; temp];
        user_points{i} = plot(temp(1), temp(2), 'o','markerfacecolor', 'b','markeredgecolor', 'b');
    end
end

if(isempty(points))
    possible_peak_locs_unsorted = locs;
else
    possible_peak_locs_unsorted = [locs; points(:,1)];
end
[possible_peak_locs, sorting_index] = sort(possible_peak_locs_unsorted);


n_new_locs = length(possible_peak_locs_unsorted);

peak_lb_ub = zeros(n_new_locs,2);


for i=1:n_new_locs
    k_loc = find(inv_locs2<possible_peak_locs_unsorted(i),1,'last');
    
    peak_lb_ub(i,:) = inv_locs2(k_loc:k_loc+1);
end

peak_lb_ub = peak_lb_ub(sorting_index,:);

end