%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm without peak analysis (= npa: np peak analysis)
% 
% 1.1 Prepare data
% 1.2 Data Preprocessing for data with sum of real part or imaginary part < 10e-6
% 2. Create predefined time constants dependent on the frequency vector
% 3. Set the regularization parameter
% 4.1 Investigating RC-elements_RC 
% 4.2 Investigating RL-elements_RC
% 4.3 Investigating R
% 4.4 Investigating L
% 4.5 Investigating C
% 5.System matrix A   
% 6.Solution vector b
% 7.1 Regularization of A 
% 7.2 Regularization of b 
% 8.1 Non-negative Least Squares Solver (NNLS)
% 8.2 Denormalize the solution
% 9. The solution of NNLS is the non-normalized h(tau)
% 10. Calculate polarization resistance ( R_pol = sum(h(tau))
% 11. Reconstruction of data
% 12. Absolute and relative error
% 13. Plot the reconstruction of data
% 14.Post-processing of the result g(tau)
% 15. Plot the whole distribution and each single distribution 
%
%
% Matthias Weiß, Yulong Zhao
%
% Last update: 20.07.2020
% last update: 29.01.2021, as function implemented, fitting functional optimization, correction of DRT
% calculation scaling problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g_DRT_RC, t, g_DRT_RL, Z_DRT, R_ohmic, C, Ax_b, x_2] ...
    = DRT_npa(Z_n, f_exp, phi, lambda)

% input:
% Z_n: input impedance (complex column vector)
% f_exp: frequency vector (column vector)
% phi: exponent for 1/(j*w*tau)^(phi) --> default = 1
% lambda: regularization parameter

% output:
% g_DRT: deconvoluted DRT
% g_DRT_peak: DRT of each peak by peak fitting
% t: discretized tau vectors = log10(tau)
% element_optimal: optimal element combination for peak fitting
% param_peak: fitting parameters for each peak
% Z_DRT: DRT reconstructed impedance
% R_ohmic: ohmic resistance of cell
% C: 1/capacitance

%% 1.1 Prepare data
Z_n = Z_n(any(Z_n, 2),:);  
f_exp = f_exp(any(f_exp, 2),:);  


re_exp = real(Z_n);                         % Extract real part of impedance
im_exp = imag(Z_n);                         % Extract imaginary part of impedance
z_exp  = re_exp + 1i*im_exp;                % Add to complex impedance

N_f = numel(f_exp);                         % Number of frequency points

f_original = f_exp;                         % Save original frequency (for plot)
z_original = z_exp;                         % Save original impedance (for plot)

%% 1.2 Data Preprocessing for data with sum of real part or imaginary part < 10e-6
% if sum(re_exp) < 10e-6 || sum(im_exp) < 10e-6
%     denorm = max([sum(re_exp),sum(im_exp)]);
%     re_exp = re_exp ./denorm;
%     im_exp = im_exp ./denorm;
% 
%     z_exp  = re_exp + 1i*im_exp;                % Add to complex impedance
% 
% else
%     denorm = 1;
% end

denorm = max(abs(z_original));
z_exp  = z_original./denorm;

%% 2. Create predefined time constants dependent on the frequency vector
N_tau = 5 * N_f;        % Number of time constants is three times number of experimental data points
fmin = min(f_exp);      % Minimum obtained experimental frequency
fmax = max(f_exp);      % Maximum obtained experimental frequency
            
tau_max = log10((1 / (fmin))) + 3;   % max(tau) +3 decades
tau_min = log10((1 / (fmax))) - 3;   % min(tau) -3 decades

tau = logspace(tau_min,tau_max,N_tau)';     % Time constants [tau_min;tau_max] 
t   = log10(tau);                           % Time constants [log10(tau_min);log10(tau_max)]

% M = length(t);
% h = 1; % mean(diff(t));
% D = ((-2).*eye(M) + diag(ones(M-1, 1), 1) + diag(ones(M-1, 1), -1))/h^2;

D = eye(length(t));
    
%% 3.1 Investigating RC-elements
omega    = 2*pi*f_exp;                              % Omega is the applied frequency
A_RC     = 1 ./ ( 1 + (1i .* omega .* tau'));       % Complex system matrix (resistive-capacitive)
Reg_RC   = lambda * D;                     % Regularization matrix (resistive-capacitive)

%% 3.2 Investigating RL-elements
omega    = 2*pi*f_exp;                              % Omega is the applied frequency
A_RL     =(1i .* omega .* tau') ./ ...
    ( 1 + (1i .* omega .* tau'));                   % Complex system matrix (resistive-inductive)
Reg_RL   = lambda * D;                     % Regularization matrix (resistive-inductive)

%% 3.3 Investigating R_ohmic
A_R      = ones(N_f,1);                             % Complex system matrix (ohmic resistance)
Reg_R    = zeros(N_tau,1);                            % Regularization matrix (ohmic resistance)

%% 3.4 Investigating L
omega    = 2*pi*f_exp;                              % Omega is the applied frequency
A_L      = 1i .* omega;                             % Complex system matrix (pure inductive)
Reg_L    = zeros(N_tau,1);                          % Regularization matrix (pure inductive)

%% 3.5 Investigating C
% C will be improved since  1/jwC is hard to model
% AC --> 1/jw * 1/C  --> algorithm is searching for 1/C
omega    = 2*pi*f_exp;                              % Omega is the applied frequency
A_C      = 1./ (1i .* omega).^phi;                       % Complex system matrix (pure capacitive)
Reg_C    = zeros(N_tau,1);                          % Regularization matrix (pure capacitive)

%% 3.System matrix A   
A_Z     = [A_RC A_RL A_R A_L A_C];  % Complex system matrix A
A_re    = real(A_Z);               	% Real part of system matrix A
A_im    = imag(A_Z);                % Imaginary part of system matrix A
A_re_im = [A_re;A_im];              % Seperated system matrix A

%% 4.Solution vector b
b_re    = real(z_exp);              % Real part of experimental data
b_im    = imag(z_exp);              % Imaginary part of experimental data
b_re_im = [b_re;b_im];              % Combination of real and imaginary part of experimental impedance

%% 5.1 Regularization of A 
Reg = [Reg_RC, Reg_RL , Reg_R, Reg_L, Reg_C];     % [Reg_RC Reg_RL  0 0 0]    
A_reg = [A_re_im; Reg];        % [A ; Reg]

%% 5.2 Regularization of b
b_reg = [b_re_im  ; zeros(N_tau,1)];       % [b ; 0]

%% 6.1 Non-negative Least Squares Solver (NNLS)
options = optimset('Display', 'off');
[x, resnorm] = lsqnonneg(A_reg,b_reg, options);   % NNLS for real and imagninare part
% x = inv(A_re_im.'*A_re_im + Reg.'*Reg)*(A_re_im.'*b_re_im);   % NNLS for real and imagninare part

%% 6.2 Denormalize the solution
x = x*denorm;

%% 7. The solution of NNLS is the non-normalized h(tau)
h_tau_RC = x(1:N_tau);                % h(tau) resistive-capacitve
h_tau_RL = x(N_tau+1:2*N_tau);        % h(tau) resistive-inductive

%% 8. Calculate polarization resistance ( R_pol = sum(h(tau))
R_pol = sum(h_tau_RC)+sum(h_tau_RL);   % The total polarization resistance RC + RL

% Normalize h(tau) by dividing by sum of h(tau)
g_tau_RC = h_tau_RC ./ R_pol;    % g(tau) = h(tau) / R_pol
g_tau_RL = h_tau_RL ./ R_pol;    % g(tau) = h(tau) / R_pol
% The integral of the distribution should lead to a 1
area_g_tau = sum(g_tau_RC) + sum(g_tau_RL);    % sum of g(tau) = 1

%% 9. Reconstruction of data
%Z_RC = (1./(1+1i*omega*tau')) * h_tau_RC;
%Z_RL = ((1i .* omega .* tau') ./( 1 + (1i .* omega .* tau')) * h_tau_RL );
Z_RC = R_pol .* ((1./(1+1i*omega*tau')) * g_tau_RC );   % R_pol * sum( g_k/(1+jw*tau_k))
Z_RL = R_pol .* ((1i .* omega .* tau') ./ ...
    ( 1 + (1i .* omega .* tau')) * g_tau_RL );
R_ohmic = x(2*N_tau + 1);
L   = x(2*N_tau + 2);
C   = x(2*N_tau + 3);
Z_L  = 1i*omega*L;
Z_C  = C./(1i*omega).^phi;

Z_DRT = Z_RC + Z_RL + R_ohmic + Z_L + Z_C;

%% 10. Absolute and relative error 
error_abs = Z_DRT-z_original;
error_rel = error_abs ./ z_original *100;
error_rel_re = real(error_rel);
error_rel_im = imag(error_rel);

%% 11. Plot the reconstruction of data

% figure('Name', 'Impedance reconstruction')
% subplot(1,2,1)
%     plot(real(z_original),-imag(z_original),'k.-.')       % Plot the experimental data
%     hold on  
%     plot(real(Z_DRT),-imag(Z_DRT),'ro')     % Plot reconstruction impedance
%     grid on
%     hold off
%         title('Electrochemical impedance spectrum','Interpreter','Latex','FontSize',14)
%         xlabel('$Z_{\mathrm{Real}}$ [$\Omega$]','Interpreter','Latex','FontSize',14)
%         ylabel('$-Z_{\mathrm{Imag}}$ [$\Omega$]','Interpreter','Latex','FontSize',14)
%         grid on  
%         legend('Experimental impedance','Reconstructed impedance','location','northwest')
% 
% subplot(1,2,2)
%     plot(log10(f_exp),error_rel_re,'-.r*')
%     hold on
%     plot(log10(f_exp),error_rel_im,':bs')
%     hold off
%         title('Relative error of reconstruction','Interpreter','Latex','FontSize',14)
%         xlabel('$\log_{10}(f_{\mathrm{exp})} [\mathrm{Hz}]$','Interpreter','Latex','FontSize',14)
%         ylabel('$\Delta_{\mathrm{rel}}$ [$\%$]','Interpreter','Latex','FontSize',14)
%         grid on  
%         legend('Relative error of the real part','Relative error of the imaginary part','location','northwest')


%% 12.Post-processing of the result g(tau) (peak analysis)

% correct DRT with regarding tau discretization
g_DRT_RC = h_tau_RC/(mean(diff(t))); % total DRT
g_DRT_RL = h_tau_RL/(mean(diff(t))); % total RL DRT

% further searching for possible hidden peaks
% diff_DRT = diff([0; h_tau_RC]);
% F_max = islocalmax(h_tau_RC);
% F_stp = islocalmin(abs(diff_DRT)); % stationary points
% 
% for i = 1:numel(F_stp)
%     if F_stp(i) == 0
%         continue
%     elseif h_tau_RC(i)  == 0 || h_tau_RC(i-1) == 0 || h_tau_RC(i+1) == 0
%         F_stp(i) = 0;
%         continue
%     elseif diff_DRT(i-1)*diff_DRT(i+1) <= 0
%         F_stp(i) = 0;
%         continue
%     elseif abs(diff_DRT(i)) >= 0.9*max(abs(diff_DRT)) 
%         F_stp(i) = 0;
%         continue
%     elseif h_tau_RC(i) <= 0.05*max(h_tau_RC)
%         F_stp(i) = 0;
%         continue
%     end
% end
% 
% % plot(diff_DRT./max(diff_DRT))
% 
% % plot previous results for assigning peak fitting elements_RC
% [h_k_RC_1,tau_k_RC_1] = findpeaks(g_DRT_RC,t, 'MinPeakHeight', 0.05*max(g_DRT_RC));    % Find peaks of the distribution
% [h_k_RL,tau_k_RL] = findpeaks(g_DRT_RL,t);    % Find peaks of the distribution
% 
% h_k_RC_2 = h_tau_RC(F_stp);
% tau_k_RC_2 = t(F_stp);
% 
% h_k_RC = [h_k_RC_1; h_k_RC_2];
% tau_k_RC = [tau_k_RC_1; tau_k_RC_2];
% 
% for i = 1:numel(h_k_RC)
%     text(tau_k_RC(i) + abs(tau_k_RC(i))*0.05, h_k_RC(i)*0.98, num2str(i))
% end
% 
% % decide whether automatic or manual or hybrid, for RL only automatic fitting is used
% prompt = 'Should automatic[A], manual[M] or hybrid[H] peak fitting be conducted?:';
% choice_peak_fitting = 'A'; % input(prompt, 's');
% 
% switch choice_peak_fitting
%     case 'A'
%         n_comb = 2^numel(h_k_RC);
%         n_comb_peak = cell(n_comb, 1);
%         for i = 1:n_comb
%             n_peak_name = dec2bin(i-1, numel(h_k_RC));
%             for j = 1:numel(h_k_RC)
%                 if n_peak_name(j) == '0'
%                     n_comb_peak{i} = [n_comb_peak{i}, 'RC,'];
%                 elseif n_peak_name(j) == '1'
%                     n_comb_peak{i} = [n_comb_peak{i}, 'RQ,'];
%                 end
%             end
%             n_comb_peak{i}(end) = '';
%         end
%         elements_RC = n_comb_peak;
%     case 'M'
%         prompt = ['assign elements_RC for each peak (#peaks = ', num2str(numel(h_k_RC)), ', separate with comma/semicolon/dot):'];
%         element = input(prompt, 's');
%         elements_RC{1} = element;
%     case 'H'
%         prompt = ['assign elements_RC for peaks (#peaks = ', num2str(numel(h_k_RC)), ', separate with comma/semicolon/dot, un-predefined element replace with ~):'];
%         element = input(prompt, 's');
%         peak_defined = strtrim(split(element, {';', '.', ','}));
%         peak_defined_index = [];
%         for n = 1:length(peak_defined)
%             if peak_defined{n} == '~'
%                 continue
%             else
%                 peak_defined_index = [peak_defined_index, n];
%             end
%         end
%         
%         n_comb = 2^numel(h_k_RC);
%         n_comb_peak = cell(n_comb, 1);
%         for i = 1:n_comb
%             n_peak_name = dec2bin(i-1, numel(h_k_RC));
%             for j = 1:numel(h_k_RC)
%                 if ismember(j, peak_defined_index)
%                     n_comb_peak{i} = [n_comb_peak{i}, [peak_defined{j}, ',']];
%                 elseif n_peak_name(j) == '0'
%                     n_comb_peak{i} = [n_comb_peak{i}, 'RC,'];
%                 elseif n_peak_name(j) == '1'
%                     n_comb_peak{i} = [n_comb_peak{i}, 'RQ,'];
%                 end
%             end
%             n_comb_peak{i}(end) = '';
%         end
%         elements_RC = unique(n_comb_peak);
%         
%     otherwise
%         disp('Either M or A must be input to choose the fitting method');
%         quit
% end
% 
% n_comb_peak = cell(numel(h_k_RL), 1);
% for i = 1:numel(h_k_RL)
%     n_peak_name = dec2bin(i-1, numel(h_k_RL));
%     for j = 1:numel(h_k_RL)
%             n_comb_peak{i} = [n_comb_peak{i}, 'RC,'];
%     end
%     n_comb_peak{i}(end) = '';
% end
% elements_RL = n_comb_peak;
% 
% [R_pol_single_RC,q_RC,gauss_RC,x_RC,element_optimal_RC] = Peak_fitting_function_DRT(g_DRT_RC,t,R_pol,h_k_RC,tau_k_RC,elements_RC);
% 
% if isempty(h_k_RL)
%     R_pol_single_RL = [];
%     q_RL = [];
%     gauss_RL = [];
%     x_RL = [];
%     element_optimal_RL = [];
% else 
%     [R_pol_single_RL,q_RL,gauss_RL,x_RL,element_optimal_RL] = Peak_fitting_function_DRT(g_DRT_RL,t,R_pol,h_k_RL,tau_k_RL,elements_RL);
% end
% 
% g_DRT_peak_RC = gauss_RC; % DRT of each peak
% g_DRT_peak_RL = gauss_RL; % DRT of each peak
% 
% %% 13. Plot the whole distribution and each single distribution
% 
% figure('Name', 'DRT deconvolution')
%     subplot(2,2,1)
%         h1 = plot(t, g_DRT_peak_RC, 'b--');
%         hold on
%         h2 = plot(t, g_DRT_RC, 'r-', 'LineWidth', 1);
%         h3 = plot(t, sum(g_DRT_peak_RC, 2), 'g-', 'LineWidth', 1);
%       	title('Deconvolution RC','Interpreter','Latex','FontSize',14)
%         xlabel('$\log_{10}(\tau)$[s]','Interpreter','Latex','FontSize',14)
%         ylabel('$q_{\log\tau}$','Interpreter','Latex','FontSize',14)
%         grid on  
%         legend([h1(1) h2 h3], 'Peak fitting', 'DRT of RC elements', 'Sum of single peaks', 'location','northwest')
%         for i = 1:numel(tau_k_RC) 
%             txt = num2str(i);
%             x_t = tau_k_RC(i);
%             y_t = h_k_RC(i) + 0.03 * max(h_k_RC);
%             text(x_t,y_t,txt)
%         end
%         hold off
%         
% if ~isempty(h_k_RL)
%     subplot(2,2,2)
%         h1 = plot(t, g_DRT_peak_RL, 'b--');
%         hold on
%         h2 = plot(t,g_DRT_RL, 'r-', 'LineWidth', 1);
%         h3 = plot(t, sum(g_DRT_peak_RL, 2), 'g-');
%        	title('Deconvolution RL','Interpreter','Latex','FontSize',14)
%         xlabel('$\log_{10}(\tau) $[s]','Interpreter','Latex','FontSize',14)
%         ylabel('$q_{\log\tau}$','Interpreter','Latex','FontSize',14)
%         grid on  
%         legend([h1(1) h2 h3], 'Peak fitting', 'DRT of RL elements', 'Sum of sinle peaks', 'location', 'northeast')
%         for i = 1:numel(tau_k_RL) 
%             txt = num2str(i+numel(tau_k_RC));
%             x_t = tau_k_RL(i);
%             y_t = h_k_RL(i) + 0.03 * max(h_k_RL);
%             text(x_t,y_t,txt)
%         end
% end
% %% 14. Error of peak fitting
% 
% total_gauss_RC = sum(g_DRT_peak_RC, 2);
% error_RC = g_DRT_RC - total_gauss_RC;
% 
% if ~isempty(h_k_RL)
%     total_gauss_RL = sum(g_DRT_peak_RL, 2);
%     error_RL = g_DRT_RL - total_gauss_RL; 
% end
% 
%     subplot(2,2,3)
%         plot(t,error_RC)
%         grid on
%         title('Error peak fitting RC','Interpreter','Latex','FontSize',14)
%         xlabel('$\log_{10}(\tau) $[s]','Interpreter','Latex','FontSize',14)
%         ylabel('$\Delta$','Interpreter','Latex','FontSize',14)
%         legend('Absolute error','location','northwest')
% 
% if ~isempty(h_k_RL)
%     subplot(2,2,4)   
%         plot(t,error_RL)
%         grid on
%         title('Error peak fitting RL','Interpreter','Latex','FontSize',14)        
%         xlabel('$\log_{10}(\tau) $[s]','Interpreter','Latex','FontSize',14)
%         ylabel('$\Delta$','Interpreter','Latex','FontSize',14)
%         legend('Absolute error','location','northeast')
% end
% %% 15. Display the results
% 
% param_peak_RC = cell(length(element_optimal_RC), 1);
% x_RC_copy = x_RC;
% for i = 1:length(element_optimal_RC)
%     switch element_optimal_RC{i}
%         case 'RC'
%             param_peak_RC{i} = x_RC_copy(1:4);
%             x_RC_copy(1:4) = [];
%         case 'RQ'
%             param_peak_RC{i} = x_RC_copy(1:3);
%             x_RC_copy(1:3) = [];
%     end
% end
% 
% disp(['R_Ohm = ', num2str(R_ohmic),' Ohm'])
% for i = 1:numel(tau_k_RC)
%     disp(['R', num2str(i) ,' = ', num2str(R_pol_single_RC(i)),' Ohm'])
% end
% 
% for i = 1:numel(tau_k_RC)
%     switch element_optimal_RC{i}
%         case 'RC'
%             disp(['Tau', num2str(i) ,' = ', num2str(10^(param_peak_RC{i}(2))),' s'])
%         case 'RQ'
%             disp(['Tau', num2str(i) ,' = ', num2str(10^(param_peak_RC{i}(2))),' s'])
%     end
% end
% 
% for i = 1:numel(tau_k_RC)
%     switch element_optimal_RC{i}
%         case 'RC'
%             disp('phi = 1')
%         case 'RQ'
%             disp(['phi', num2str(i) ,' = ', num2str(param_peak_RC{i}(3))])
%     end
% end
% 
% 
% param_peak_RL = cell(length(element_optimal_RL), 1);
% if ~isempty(h_k_RL)
%     x_RL_copy = x_RL;
%     for i = 1:length(element_optimal_RL)
%             param_peak_RL{i} = x_RL_copy(1:4);
%             x_RL_copy(1:4) = [];
%     end
% 
%     for i = 1:numel(tau_k_RL)
%         disp(['R_RL', num2str(i) ,' = ', num2str(R_pol_single_RL(i)), ' Ohm'])
%         disp(['tau_RL', num2str(i), ' = ', num2str(10^(param_peak_RL{i}(2))), ' s'])
%     end
%     
% end
end

