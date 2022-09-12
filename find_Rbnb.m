function [cost, R, y_hat, res, R_min, res_der, Rmin_indices, int1, cost_in_vec] = find_Rbnb(A,y, is_lasso, cand_info, disable_vec)

persistent myopt F max_y

if(isempty(myopt) || isempty(max_y))
    myopt =  optimset('display','off'); %,'Algorithm', 'active-set'
    max_y = max(y);
    
    
    %     sdpv_x = sdpvar(n2,1);
    %     sdpv_A = sdpvar(n1,n2,'full');
    %
    %     const = [0 <= sdpv_x <= 1.2*max_y];
    %
    %     obj = norm(sdpv_A*sdpv_x - y,2)^2;
    %
    %     F = optimizer(const, obj, sdpsettings('sovler','gurobi'), sdpv_A, sdpv_x);
    %
    %
    
    
end


if(nargin<3)
    is_lasso = false;%true;% false
end

if(nargin<4)
    is_same_mu = false;
else
    is_same_mu = cand_info.is_same_mu;
end

if(nargin<5)
    disable_vec_flat = false(size(A,2),1);
else
    disable_vec_flat = R_flatten(disable_vec,cand_info);
end


A(:,any(isnan(A))) = 0;
A(:,disable_vec_flat') = 0;

lb = zeros(size(A,2),1);
ub = 1.2*max_y*ones(size(A,2),1);




if(is_lasso)
    R = yalmip_intquad(A, y, lb, ub);
else
   % R = lsqlin(A,y,[],[],[],[], lb, ub, [], myopt);
   R = lsqnonneg(A,y,myopt);
end

if(is_same_mu)
    %cost_double = 0;
    R_shaped = R_shaper(R, cand_info);

    %    R_shaped = reshape(R,[],cand_info.n_cand_eff);
    [R_min, Rmin_indices] = mink(R_shaped,size(R_shaped,1)-1,1);%min(R_shaped,[],1);
    cost_double = sum(R_min(:));
    
    % ---- for ignoring small ones ----- 
    R_shaped2 = R_shaped;
    for ix = 1:size(Rmin_indices,2)
        R_shaped2(Rmin_indices(:,ix),ix) = 0;
    end
    
    
    R2 = R_flatten(R_shaped2, cand_info);
    
    % ------------------------------------
    
else
    cost_double = 0;
end

%R =   lsqlin(A,y,[],[],[],[], lb, ub, [], myopt); % quadprog_lasso(A, y, lb, ub); %2x slower than lsqlin but much better than sdpvar.
%lsqlin(A,y,[],[],[],[],zeros(size(A,2),1), 1.2*max_y*ones(size(A,2),1), [], myopt);
%toc
%R = A\y;
y_hat = A*R2;
y_hat_all = A.*R2';

int1 = trapz(y_hat_all)/size(A,1);
int1(R2'<=0.01) = 0;



% plot(y); hold on; plot(y_hat,'--');

res = y - y_hat;
%res_der = diff(y) - diff(y_hat); 
res_der = 0;
%cost_of_max = sum(mink(R,max(length(R)-9,0))); % put a cost more than 9 R values.

cost_in_vec = [10*res/length(res); sqrt(R(:))/length(R); 5*int1(:)]; 



%cost1 = norm(res,2)^2/length(res);
%cost2 = 0.01*norm(R,2)^2/length(R) + 500*cost_double/cand_info.n_cand;%norm(R,1); % + 200*cost_double/cand_info.n_cand;% + 50*cost_of_max;
%cost3 = (norm(res_der,2)/length(y))^2;

%cost_area = norm(int1,2)^2;
cost  = norm(cost_in_vec,2)^2;% 100*cost1  + 0.01*cost2 + cost_area;% + 0.01*cost3;%+ 0.2*max(res)/max_y + 0.5*norm(diff(y) - diff(y_hat),2)^2/max(abs(diff(y)));

end