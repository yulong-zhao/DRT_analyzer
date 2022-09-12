function [last_param, ys] = printModelparam(x, cand, param, cand_info, y, y_orj)
trigo_kernel_old = @(x, phi, mu)  sin(phi*pi)./((cosh(phi.*(x-mu) ) +cos(phi*pi)  ) );
trigo_kernel_old_fun  = @(x,p) trigo_kernel_old(x, p(1,:), p(2,:));
    
is_same_mu = cand_info.is_same_mu;
n_space = cand_info.n_space;

A = zeros(length(x), n_space);
A2 = zeros(length(x), n_space);

last_param = {};

bnd2 = 0;
param_count = 0;
if(is_same_mu)
    mu_all = reshape(param(1:cand_info.n_peaks),1,[]);
    param_count = cand_info.n_peaks;
end

bnds = [];

A_other = 0;
for i=1:cand_info.n_cand
    
    bnd1 = bnd2 + 1;
    bnd2 = bnd1 - 1 + cand(i).n_funs;
    
    bnds = [bnds; bnd1, bnd2];
    
    this_params = param(param_count+1 : param_count+cand(i).n_p_tot);
    
    
    
    reshaped_param = reshape(this_params, cand(i).n_param, []);
    
    %
    reshaped_param = [reshaped_param; mu_all(cand_info.indices{i})];
    last_param{end+1} = reshaped_param;
    if(~isempty(reshaped_param))
        if i==2 % Use old trigo kernel to get actual R.
            A(:, bnd1:bnd2)= cand(i).fun(x, reshaped_param);
            A_other = trigo_kernel_old_fun(x, reshaped_param);
            A2(:, bnd1:bnd2)= A_other;
        else
            A(:, bnd1:bnd2)= cand(i).fun(x, reshaped_param); %trigo_kernel_old_fun(x, reshaped_param);
            A2(:, bnd1:bnd2)= cand(i).fun(x, reshaped_param);
        end
        
    end
    
    param_count = param_count+cand(i).n_p_tot;
end
threshold =0.01;
[cost, R_opt,  y_hat, res] = find_R(A, y, 0, cand_info);


max_trigo = max(A_other);

R_opt_clean = R_opt;

R_opt_clean(R_opt<threshold) = 0;

R_opt_adjusted = R_opt_clean;

param_names = { {'A  ', 'R  ', 'W  ', 'mu '}, {'A  ', 'R  ', 'phi', 'mu '}, {'A  ', 'R  ', 's  ', 'W  ', 'mu '}};

for i=1:cand_info.n_cand
    R_new = vec(R_opt_clean(bnds(i,1):bnds(i,2)))'*max(y_orj);
    
    if(i==2)
        R_new = R_new./max_trigo;
    end
    
    R_opt_adjusted(bnds(i,1):bnds(i,2)) = R_new;
    
    last_param{i} = [R_new  ;last_param{i}];
end


ys = A2.*R_opt_adjusted';
areas = trapz(x,ys,1);

for i=1:cand_info.n_cand
    last_param{i} = [areas(bnds(i,1):bnds(i,2)); last_param{i}];
end


for i=1:cand_info.n_cand
    fprintf('Candidate : %d\n',i);
    for j=1:size(last_param{i},1)
        if(j>2)
            fmt=[param_names{i}{j}, '=' repmat(' %4.6f\t\t',1,numel(last_param{i}(j,:))) '\n'];
        else
            fmt=[param_names{i}{j}, '=' repmat(' %4.3e\t\t',1,numel(last_param{i}(j,:))) '\n'];
        end
        
        fprintf(fmt,last_param{i}(j,:));
        
    end
    fprintf('\n');
end


% figure; 
% title('Original data');
% plot(x, y_orj); hold on;
% plot(x, A2*R_opt_adjusted, '--');
% plot(x, ys);
% legend('orj', 'estimated');

end