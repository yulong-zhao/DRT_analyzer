
if( n_funs > n_locs_try )
    n_diff = (n_funs - n_locs_try);
    warm.locs =  round(linspace(1,length(x_try),n_funs)); % locs; % %round(logspace(0,log10(length(x)),n_funs));%
    warm.x_loc = x_try(warm.locs);
    warm.W = w_try; %rand(size(warm.x_loc));
    warm.y_loc = y_try(warm.locs);
    warm.s = zeros(1,n_funs);
    
    i_cand = 1;
    
    warm.bnds_low{i_cand} = [ cand(i_cand).warm.low.W ; cand(i_cand).warm.low.mu];
    warm.bnds_low{i_cand} = warm.bnds_low{i_cand}(:);
    
    warm.bnds_high{i_cand} = [ cand(i_cand).warm.high.W ; cand(i_cand).warm.high.mu];
    warm.bnds_high{i_cand} = warm.bnds_high{i_cand}(:);
    
    temp = [[2*warm.W(:)', rand(1, n_diff)]; warm.x_loc(:)'];
    warm.alls{i_cand} = temp(:);
    
    
    
    i_cand = 2;
    
    warm.bnds_low{i_cand} = [ cand(i_cand).warm.low.phi ; cand(i_cand).warm.low.mu];
    warm.bnds_low{i_cand} = warm.bnds_low{i_cand}(:);
    
    warm.bnds_high{i_cand} = [ cand(i_cand).warm.high.phi ; cand(i_cand).warm.high.mu];
    warm.bnds_high{i_cand} = warm.bnds_high{i_cand}(:);
    
    temp = [[2*warm.W(:)', rand(1, n_diff)]; warm.x_loc(:)']; %[warm.y_loc(:)';  warm.W(:)'; warm.x_loc(:)'];
    
    warm.alls{i_cand} = temp(:);
    
    
    i_cand = 3;
    
    warm.bnds_low{i_cand} = [ cand(i_cand).warm.low.S; cand(i_cand).warm.low.W ; cand(i_cand).warm.low.mu];
    warm.bnds_low{i_cand} = warm.bnds_low{i_cand}(:);
    
    warm.bnds_high{i_cand} = [cand(i_cand).warm.high.S; cand(i_cand).warm.high.W ; cand(i_cand).warm.high.mu];
    warm.bnds_high{i_cand} = warm.bnds_high{i_cand}(:);
    
    temp = [warm.s; [2*warm.W(:)', rand(1, n_diff)]; warm.x_loc(:)']; %[warm.y_loc(:)';  warm.W(:)'; warm.x_loc(:)'];
    
    warm.alls{i_cand} = temp(:);
    
    
    
    warm.all = vertcat(warm.alls{which_funs}); %repmat(vec([warm.y_loc(:)';  warm.W(:)'; warm.x_loc(:)']),n_cand,1);
    
    warm.bnds_low_all  =  vertcat(warm.bnds_low{which_funs});
    warm.bnds_high_all =  vertcat(warm.bnds_high{which_funs});
    
    
    
    PSO_InitialSwarmMatrix  = unifrnd(repmat(warm.bnds_low_all',N_PSO_mat,1), repmat(warm.bnds_high_all',N_PSO_mat,1));
    
    warm.bnds_low_all_ga  = [warm.bnds_low_all; zeros(n_funs,1)];
    warm.bnds_high_all_ga = [warm.bnds_high_all; n_cand*ones(n_funs,1)];
    
    GA_InitialSwarmMatrix  = unifrnd(repmat(warm.bnds_low_all_ga',N_PSO_mat,1), repmat(warm.bnds_high_all_ga',N_PSO_mat,1));
    
    GA_InitialSwarmMatrix(:,(n_param+1):n_ga) = round(GA_InitialSwarmMatrix(:,(n_param+1):n_ga));
    
else
    warm.locs =  round(linspace(1,length(x),n_funs)); % locs; % %round(logspace(0,log10(length(x)),n_funs));%
    warm.x_loc = x(warm.locs);
    warm.W = w_try; %rand(size(warm.x_loc));
    warm.y_loc = y(warm.locs);
    
    temp = [2*warm.W(:)'; warm.x_loc(:)'];
    warm.alls{1} = temp(:);
    % temp = [warm.W(:)'; warm.x_loc(:)']; %[warm.y_loc(:)';  warm.W(:)'; warm.x_loc(:)'];
    % warm.alls{1} = temp(:);
    warm.all = vertcat(warm.alls{which_funs}); %repmat(vec([warm.y_loc(:)';  warm.W(:)'; warm.x_loc(:)']),n_cand,1);
    
end