function configs = sim_config(state, sig_q)
configs.state = state;
configs.P = diag([25,25,25,5e-6,5e-6,5e-6]);
% configs.Q = (sig_w^2)*...
%     [((dt^3)/3), 0, 0, ((dt^2)/2), 0, 0;
%     0, ((dt^3)/3), 0, 0, ((dt^2)/2), 0;
%     0, 0, ((dt^3)/3), 0, 0, ((dt^2)/2);
%     ((dt^2)/2), 0, 0, dt, 0, 0;
%     0, ((dt^2)/2), 0, 0, dt, 0;
%     0, 0, ((dt^2)/2), 0, 0, dt];

configs.sig_q = sig_q;

% configs.Q = zeros(6,6);

configs.R = diag([(10e-3)^2, (1e-3)^2, (1e-3)^2, (1e-3)^2]);

end