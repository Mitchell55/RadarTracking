classdef BatchFilter < handle
    properties
        R
        P
    end
    
    methods
        % add a constructor
        function obj = BatchFilter(configs)
            obj.R = configs.R;
        end
        
        function [initial_state, filter_residuals, eciMeas] = process_meas(obj, meas, t, lla, vel)
            % convert to eci
            eciMeas = obj.convertToEci(meas, t, lla);
            initial_state = obj.getInitialState(eciMeas, t, vel);
            filter_residuals = obj.getFilterResiduals(initial_state, eciMeas, t, vel);
        end
        
        function eciMeas = convertToEci(obj, meas, t, lla)
            eciMeas = nan(3, size(meas, 2));
            for i = 1:size(meas,2)
                [eciMeas(:,i), ~ ] = gnc.getstate(meas(:,i), lla(i,:), t(i));
            end

        end

        function initial_state = getInitialState(obj, eciMeas, t, vel)
            % First we need to stack all the stat transition matricies
            options = odeset('RelTol',1e-13,'AbsTol',1e-15);
            % F_tot = eye(6);
            h = [1, 0, 0, 0, 0, 0;...
                  0, 1, 0, 0, 0, 0;...
                  0, 0, 1, 0, 0, 0];
            H_tot = h;
            Y_tot = eciMeas(:,1);
            R_total = obj.R;
            L = h*eye(6);
            for i = 2:numel(t)
                % prop forward to get the velocity
                dt = seconds( t(i) - t(1));
                Phi_flat = reshape(eye(6),6^2,1);
                state = [eciMeas(:,1); vel; Phi_flat];
                [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
                    [0, dt],state,...
                    options);
                current_vel = final_state(end,4:6)';
                dt = seconds(t(1) - t(i));
                state = [eciMeas(:,i); current_vel; Phi_flat];
                % make the measurements go back in time
                [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
                    [0, dt],state,...
                    options);
                % concat the backward in time measurements
                Y_tot = [Y_tot; final_state(end,1:3)'];
                % concat the meansurement matricies
                H_tot = [H_tot; h];
                % H_tot = blkdiag(H_tot, h);
                % concat the L matrix
                % has to be inverse because we went back in time.
                L = [L; h* inv(reshape(final_state(end,7:end), [6, 6])) ];
                % concat the measurements
                % Y_tot = [Y_tot; eciMeas(:,i)];
                % concat the R matricies
                % R_total = [R_total; obj.R];
                R_total = blkdiag(R_total, obj.R);
            end
            % solve for the initial state
            initial_state = inv( (H_tot)' * inv(R_total) *  (H_tot) ) * (H_tot)' * inv(R_total) * Y_tot;
            obj.P = inv(H_tot' * inv(R_total) *H_tot);

        end
        
        
        function filter_residuals = getFilterResiduals(obj, initial_state, eciMeas, t, vel)
            filter_residuals = initial_state(1:3,1) - eciMeas(:,1);
            options = odeset('RelTol',1e-13,'AbsTol',1e-15);
            for i = 2: numel(t)
                dt = seconds( t(i) - t(1) ) ;
                Phi_flat = reshape(eye(6),6^2,1);
                state = [initial_state(1:3,1); vel;  Phi_flat];
                [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
                    [0, dt],state,...
                    options);
                filter_residuals = [filter_residuals, final_state(end,1:3)' - eciMeas(:, i)];

            end


        end
    
    
    end
    

    
    
    
    
end