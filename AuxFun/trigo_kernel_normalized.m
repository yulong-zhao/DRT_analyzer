function norm_trig_out = trigo_kernel_normalized(x, phi, mu)
norm_trig_out = sin(phi*pi)./((cosh(phi.*(x-mu) ) +cos(phi*pi)  ) );

norm_trig_out = norm_trig_out./max(norm_trig_out); % So the peak is always one. 

end