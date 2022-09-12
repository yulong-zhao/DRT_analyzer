function [fitresult, gof] = createWarmStart(x_try, y_try, cand, cand_info, opt, warm)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 24-Jul-2021 04:26:49


wcr.shaped_x0 = zeros(3,cand_info.n_peaks);
wcr.shaped_x0(2,:) = warm.all_x0(1:cand_info.n_peaks);
wcr.shaped_x0(3,:) = mean(cand(1).n_limits_warm(1,:))/2;
wcr.shaped_x0(1,:) = interp1(x_try,y_try,wcr.shaped_x0(2,:));

wcr.shaped_lb = zeros(3,cand_info.n_peaks);
wcr.shaped_lb(1,:) = 0;
wcr.shaped_lb(2,:) = opt.lb(1:cand_info.n_peaks);
wcr.shaped_lb(3,:) = cand(1).n_limits(1,1);

wcr.shaped_ub = zeros(3,cand_info.n_peaks);
wcr.shaped_ub(1,:) = 1;
wcr.shaped_ub(2,:) = opt.ub(1:cand_info.n_peaks);
wcr.shaped_ub(3,:) = cand(1).n_limits(1,2);



%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x_try, y_try );

% Set up fittype and options.
ft = fittype( 'gauss7' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMaxChange = 0.01;
opts.Display = 'Off';
opts.Lower = wcr.shaped_lb(:)'; % R, mu, sigma
opts.Upper = wcr.shaped_ub(:)'; %[1 Inf Inf 1 Inf Inf 1 Inf Inf 1 Inf Inf 1 Inf Inf 1 Inf Inf 1 Inf Inf];
opts.MaxFunEvals = 5000;
opts.MaxIter = 5000;
opts.Normalize = 'on';
opts.StartPoint =  wcr.shaped_x0(:)';

%[0.0516229108639802 0.167530859202258 0.0298326531523249 0.0495783501778315 0.264806841964859 0.0334026240555779 0.0448278299202675 0.113488646556368 0.0347874313354089 0.0396501232266082 0.21076462931897 0.0421562090569843 0.0384418969481636 0.318849054610749 0.0508715844223953 0.0349653895793103 0.567443232781841 0.0716463719149355 0.0288850141260961 0.0594464339104789 0.0530445570902786];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on

end







