InitialPoints.n = 1000; 
InitialPoints.loc_perturb = 0.05;

InitialPoints.loc_lb = warm.all_x0(1:cand_info.n_peaks) - InitialPoints.loc_perturb;
InitialPoints.loc_ub = warm.all_x0(1:cand_info.n_peaks) + InitialPoints.loc_perturb;

InitialPoints.lb = [InitialPoints.loc_lb; opt.lb(cand_info.n_peaks+1:end)];
InitialPoints.ub = [InitialPoints.loc_ub; opt.ub(cand_info.n_peaks+1:end)];    



InitialPoints.all = [warm.all_x0'; unifrnd(repmat(InitialPoints.lb',InitialPoints.n-1,1), repmat(InitialPoints.ub',InitialPoints.n-1,1))];

InitialPoints.all = max(InitialPoints.all, opt.lb');
InitialPoints.all = min(InitialPoints.all, opt.ub');


InitialPoints.points = CustomStartPointSet(InitialPoints.all);

