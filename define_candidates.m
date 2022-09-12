gauss_kernel = @(x, W,   mu)  exp(- ( (x-mu)./W ).^2);
trigo_kernel_old = @(x, phi, mu)  sin(phi*pi)./((cosh(phi.*(x-mu) ) +cos(phi*pi)  ) );
trigo_kernel = @(x, phi, mu)  (1+cos(phi*pi))./((cosh(phi.*(x-mu) ) +cos(phi*pi)  ) );                               %@trigo_kernel_normalized;
nsym_gauss_k_wo_sign = @(x, s, sigma, mu) exp(- ((x-mu).^2.*(1+s.^2) + 2*sqrt((x-mu).^2).*(x-mu).*s)./(2*sigma.^2));

nsym_gauss_k = @(x, s, sigma, mu) nsym_gauss_k_wo_sign(x, s, sigma, mu); %exp(- ((x-mu).*(1 + sign(x-mu).*s)).^2./(2*sigma.^2));  % -2<=s<=2


sensitivity.dx = 1;
sensitivity.dw = 1;


cand_info.n_cand_eff = nnz(cand_info.n_funs_all);

i_cand = 1;
cand(i_cand).fun            = @(x,p) gauss_kernel(x, p(1,:), p(2,:));
cand(i_cand).n_limits       = [0.13, 0.92; min(x_try), max(x_try)]; % Important. 
cand(i_cand).n_limits_warm  = [0.3, 0.5; min(x_try), max(x_try)];


% %
i_cand = 2;
cand(i_cand).fun            = @(x,p) trigo_kernel(x, p(1,:), p(2,:));
cand(i_cand).n_limits       = [0.81, 0.98; min(x_try), max(x_try)];
cand(i_cand).n_limits_warm  = [0.85, 0.96; min(x_try), max(x_try)];

%
i_cand = 3;
cand(i_cand).fun            = @(x,p) nsym_gauss_k(x, p(1,:), p(2,:), p(3,:));
cand(i_cand).n_limits       = [-5, 5; 0.09, 0.8; min(x_try), max(x_try)];
cand(i_cand).n_limits_warm  = [-0.1, 0.1; 0.2, 0.4; min(x_try), max(x_try)];



for i_cand = 1:cand_info.n_cand
    cand(i_cand).n_param        = get_N_param(cand(i_cand).fun);
    cand(i_cand).n_funs         = cand_info.n_funs_all(i_cand);
    cand(i_cand).n_diff         = 0;%cand(i_cand).n_funs - n_locs_try;
    cand(i_cand).n_p_tot        = cand(i_cand).n_param*cand(i_cand).n_funs;
    
    
    cand(i_cand).R_limits       = [0, 1.2];
    
%     if(i_cand==3)
%         cand(i_cand).warm.x0        = [[zeros(size(w_try(:)')); w_try(:)'; locs_try(:)'], ...
%             unifrnd( repmat(cand(i_cand).n_limits_warm(:,1), 1, cand(i_cand).n_diff), ...
%             repmat(cand(i_cand).n_limits_warm(:,2), 1, cand(i_cand).n_diff))];
%     else
% 
%         cand(i_cand).warm.x0        = [[w_try(:)'; locs_try(:)'], ...
%             unifrnd( repmat(cand(i_cand).n_limits_warm(:,1), 1, cand(i_cand).n_diff), ...
%             repmat(cand(i_cand).n_limits_warm(:,2), 1, cand(i_cand).n_diff)  )];
%         
%     end
    
    cand(i_cand).n_limits_eff = cand(i_cand).n_limits;
    cand(i_cand).n_limits_warm_eff = cand(i_cand).n_limits_warm;
  %  cand(i_cand).warm.x0_eff  = cand(i_cand).warm.x0;
    
    if(cand_info.is_same_mu) 
        % if it is same mu then remove the mu parameters for the functions 
        cand(i_cand).n_limits_eff(end,:) = [];
        cand(i_cand).n_limits_warm_eff(end,:) = [];
    %    cand(i_cand).warm.x0_eff(end,:)  = [];
        cand(i_cand).n_param         = cand(i_cand).n_param - 1;
        cand(i_cand).n_p_tot         = cand(i_cand).n_param*cand(i_cand).n_funs;
    end
        
        
    
    
    
    %    cand(i_cand).warm.low.W     = [2*w_try(:)' - sensitivity.dw, repmat(cand(i_cand).n_limits(1,1), 1, cand(i_cand).n_diff)];
    %    cand(i_cand).warm.high.W    = [2*w_try(:)' + sensitivity.dw, repmat(cand(i_cand).n_limits(1,2), 1, cand(i_cand).n_diff)];
    %    cand(i_cand).warm.low.mu    = [locs_try'   - sensitivity.dx, repmat(cand(i_cand).n_limits(2,1), 1, cand(i_cand).n_diff)];
    %    cand(i_cand).warm.high.mu   = [locs_try'   + sensitivity.dx, repmat(cand(i_cand).n_limits(2,2), 1, cand(i_cand).n_diff)];
end



cand_info.n_param_per_cand = [cand(:).n_p_tot];
cand_info.n_param = sum(cand_info.n_param_per_cand);

if(cand_info.is_same_mu)
    cand_info.n_param = cand_info.n_param + cand_info.n_peaks;
end

cand_info.n_space = sum([cand(:).n_funs]);



% cand(i_cand).warm.x0        = [[w_try(:)'; locs_try(:)'], ...
%                               unifrnd( repmat(cand(i_cand).n_limits_warm(:,1), 1, cand(i_cand).n_diff), ...
%                                        repmat(cand(i_cand).n_limits_warm(:,2), 1, cand(i_cand).n_diff))];
%
%
% cand(i_cand).warm.low.phi   = [2*w_try(:)' - sensitivity.dw, repmat(cand(i_cand).n_limits(1,1), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.high.phi  = [2*w_try(:)' + sensitivity.dw, repmat(cand(i_cand).n_limits(1,2), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.low.mu    = [locs_try'   - sensitivity.dx, repmat(cand(i_cand).n_limits(2,1), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.high.mu   = [locs_try'   + sensitivity.dx, repmat(cand(i_cand).n_limits(2,2), 1, cand(i_cand).n_diff)];
%
% cand(i_cand).warm.low.phi     = max(cand(i_cand).warm.low.phi,  cand(i_cand).n_limits(1,1));
% cand(i_cand).warm.high.phi    = min(cand(i_cand).warm.high.phi, cand(i_cand).n_limits(1,2));
%
% cand(i_cand).warm.low.mu     = max(cand(i_cand).warm.low.mu,  cand(i_cand).n_limits(2,1));
% cand(i_cand).warm.high.mu    = min(cand(i_cand).warm.high.mu, cand(i_cand).n_limits(2,2));
%
%
%

%
%
% cand(i_cand).warm.low.S     = repmat(cand(i_cand).n_limits(1,1), 1, cand(i_cand).n_funs);
% cand(i_cand).warm.high.S    = repmat(cand(i_cand).n_limits(1,2), 1, cand(i_cand).n_funs);
% cand(i_cand).warm.low.W     = [2*w_try(:)' - sensitivity.dw, repmat(cand(i_cand).n_limits(2,1), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.high.W    = [2*w_try(:)' + sensitivity.dw, repmat(cand(i_cand).n_limits(2,2), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.low.mu    = [locs_try'   - sensitivity.dx, repmat(cand(i_cand).n_limits(3,1), 1, cand(i_cand).n_diff)];
% cand(i_cand).warm.high.mu   = [locs_try'   + sensitivity.dx, repmat(cand(i_cand).n_limits(3,2), 1, cand(i_cand).n_diff)];

% cand(i_cand).warm.low.W     = max(cand(i_cand).warm.low.W,  cand(i_cand).n_limits(2,1));
% cand(i_cand).warm.high.W    = min(cand(i_cand).warm.high.W, cand(i_cand).n_limits(2,2));


% if(is_same_mu)
%     cand(i_cand).warm.low.mu     = [];% max(cand(i_cand).warm.low.mu,  cand(i_cand).n_limits(3,1));
%     cand(i_cand).warm.high.mu    = [];% min(cand(i_cand).warm.high.mu, cand(i_cand).n_limits(3,2));
%     cand(i_cand).n_param         = cand(i_cand).n_param - 1;
%     cand(i_cand).n_p_tot         = cand(i_cand).n_param*cand(i_cand).n_funs;
%     cand(i_cand).n_limits(end,:) = [];
%     cand(i_cand).n_limits_warm(end,:) = [];
%     cand(i_cand).warm.x0(end,:) = [];
% else
%     cand(i_cand).warm.low.mu     = max(cand(i_cand).warm.low.mu,  cand(i_cand).n_limits(3,1));
%     cand(i_cand).warm.high.mu    = min(cand(i_cand).warm.high.mu, cand(i_cand).n_limits(3,2));
% end
