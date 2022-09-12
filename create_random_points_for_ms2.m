InitialPoints.n = 5000; 
InitialPoints.loc_perturb = 0.5;

InitialPoints.loc_lb = warm.all_x0(1:cand_info.n_peaks) - InitialPoints.loc_perturb;
InitialPoints.loc_ub = warm.all_x0(1:cand_info.n_peaks) + InitialPoints.loc_perturb;

InitialPoints.lb = [InitialPoints.loc_lb; opt.lb(cand_info.n_peaks+1:end)];
InitialPoints.ub = [InitialPoints.loc_ub; opt.ub(cand_info.n_peaks+1:end)];    

InitialPoints.withconstantloc = [repmat(warm.all_x0(1:cand_info.n_peaks)',InitialPoints.n-1,1), unifrnd(repmat(opt.lb(cand_info.n_peaks+1:end)',InitialPoints.n-1,1), repmat(opt.ub(cand_info.n_peaks+1:end)',InitialPoints.n-1,1))];
InitialPoints.withperturbedloc = unifrnd(repmat(InitialPoints.lb',InitialPoints.n-1,1), repmat(InitialPoints.ub',InitialPoints.n-1,1));

%InitialPoints.all = [warm.all_x0'; InitialPoints.withperturbedloc];

InitialPoints.all = [warm.all_x0'; InitialPoints.withconstantloc];

InitialPoints.all = max(InitialPoints.all, opt.lb');
InitialPoints.all = min(InitialPoints.all, opt.ub');



% test_x0 = [-4.2829; -3.4103; -2.8140; -2.1322; -0.7053; 0.1145; 1.2486;...
%          -0.4114; 0.1275; ...
%          0.0423; 0.2994; ...
%          -0.0574; 0.3059; ...
%          -0.0562; 0.3568; ...
%          0.4214; 0.1831; ...
%          -0.2271; 0.2487; ...
%          -0.1777; 0.3002];
% 
% 
% InitialPoints.all = [test_x0'; InitialPoints.all];

     




InitialPoints.points = CustomStartPointSet(InitialPoints.all);




     

