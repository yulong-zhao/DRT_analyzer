%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the optimal lambda using discrepancy method

function lambda_opt = discrepancy_lambda(z_exp, f_exp)

N_lambda = 100;
lambda_set = logspace(-3, 1, N_lambda); % 

parfor i = 1:N_lambda
    [x_re(:, i), x_im(:, i)] ...
        = DRT_npa_lambda_opt(z_exp, f_exp, 1, lambda_set(i));
    delta_re_im(i) = norm(x_re(:, i) - x_im(:, i));
end

[~, idx_opt] = min(delta_re_im);
lambda_opt = lambda_set(idx_opt);

end