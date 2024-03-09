function configs = real_data_config(state, sig_q)
configs.state = state;
configs.P = diag([25,25,25,5e-6,5e-6,5e-6]);

configs.sig_q = sig_q;

% configs.Q = zeros(6,6);

configs.R = diag([(10e-3)^2, (1e-3)^2, (5e-3)^2, (5e-3)^2]);

end