function [Z_peak, fit_param] = CPE_calculate(f_exp, t, g_peak)

omega = 2*pi.*f_exp;
tau = 10.^t';


Z_peak     = trapz(t', g_peak'./(1+1j.*omega*tau), 2);
f_obj_peak = @(x) sum(abs((real(Z_peak) - x(1)).^2 + (imag(Z_peak) - x(2)).^2 - (x(3))^2))*1e3;
x0_peak    = [(max(real(Z_peak)) + min(real(Z_peak)))/2, -1e-3, (max(real(Z_peak)) - min(real(Z_peak)))/2];
options    = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, 'FunctionTolerance', 1e-6, ...
    'StepTolerance', 1e-6, 'Display', 'none');  
[x_RC_peak, ~] = lsqnonlin(f_obj_peak, x0_peak, [], [], options);

angle = asin(x_RC_peak(2)/x_RC_peak(3)) + pi/2;
phi_CPE = angle/(pi/2);
cx = x_RC_peak(1);
cy = x_RC_peak(2);
r  = x_RC_peak(3);

fit_param = [cx, cy, r, phi_CPE];

end