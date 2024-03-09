classdef EKF < handle
    properties
        state
        P
        sig_q
        R
        current_time
    end
    
    methods
        % add a constructor
        function obj = EKF(configs, current_time)
            obj.state = configs.state;
            obj.P = configs.P;
            obj.sig_q = configs.sig_q;
            obj.R = configs.R;
            obj.current_time = current_time;
        end
        
        
        function filter_resid = process_meas(obj, meas, t, lla)
            obj.predict(t);
            obj.current_time = t;
            filter_resid = obj.meas_update(meas, lla, t);
            
        end
        
        
        function predict(obj, t_prop)
            % need the F matrix
            % F = obj.getF(obj.state);
            
            % ekf uses the nonlinear propagation
            options = odeset('RelTol',1e-13,'AbsTol',1e-15);
            dt = seconds(t_prop - obj.current_time);
            if dt < 0
                disp('Error: went back in time');
            end
            Phi_flat = reshape(eye(6),6^2,1);
            state = [obj.state;Phi_flat];
            [prop_times,final_state] = ode45(@(t,final_state) gnc.state_dyn(t,final_state),...
                [0, dt],state,...
                options);
            F = reshape(final_state(end, 7:end)', 6,6 ) ;
            obj.state = final_state(end,1:6)';

            Q = gnc.getQ(obj.sig_q, dt);
            
            obj.P = F * obj.P * F' + Q;
          
        end
        
        function filter_resid = meas_update(obj, meas, lla, time)
            % get H tilde
            H = obj.getH(obj.state, lla, time);
            % predected observation 
            predicted_radar_observation = gnc.gen_observation_fn(obj.state, lla, time );
            % range, range rate, az (deg), el (deg)
            ey = meas - predicted_radar_observation;
            % Calculate the Kalman Gain
            k = obj.P*H' * inv(H * obj.P * H' + obj.R);
            % state update
            obj.state = obj.state + k*ey;
            % covariance update - P_temp = (eye(size(P))-K*H_ekf)*P*(eye(size(P))-K*H_ekf)' + K*model.R*K'; % Joseph-Bucy form!!!
            obj.P = (eye(6) - k*H)* obj.P * (eye(6) - k*H)' + k*obj.R*k';
            filter_resid = ey;
        end
        
        function F = getF(obj,state)
            % set up constants
            x = state(1);
            y = state(2);
            z = state(3);
            mu = 398600.435507;
            J2 = 1.75553e10;
            rE = 6378.1;
            
            
            F = [0,0,0, 1, 0, 0;...
                 0,0,0, 0, 1, 0;...
                 0,0,0,0,0,1;...
                 (mu*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) + (mu*((3*J2*rE^2*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2 - (24*J2*rE^2*x^2*z^2)/(x^2 + y^2 + z^2)^4 - (4*J2*rE^2*x^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (2*mu*x*((3*J2*rE^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*x*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x^2*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), - (mu*((24*J2*rE^2*x*y*z^2)/(x^2 + y^2 + z^2)^4 + (4*J2*rE^2*x*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*y*((3*J2*rE^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*x*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*x*((3*J2*rE^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x*y*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), (mu*x*((J2*rE^2*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (J2*rE^2*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*z*((3*J2*rE^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*x*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*((18*J2*rE^2*x*z^3)/(x^2 + y^2 + z^2)^4 - (J2*rE^2*x*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^2 - (6*J2*rE^2*x*z)/(x^2 + y^2 + z^2)^3 + (4*J2*rE^2*x*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (3*mu*x*z*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0;...
                 - (mu*((24*J2*rE^2*x*y*z^2)/(x^2 + y^2 + z^2)^4 + (4*J2*rE^2*x*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*y*((3*J2*rE^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*x*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*x*((3*J2*rE^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x*y*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), (mu*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) + (mu*((3*J2*rE^2*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2 - (24*J2*rE^2*y^2*z^2)/(x^2 + y^2 + z^2)^4 - (4*J2*rE^2*y^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (2*mu*y*((3*J2*rE^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*y^2*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), (mu*y*((J2*rE^2*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (J2*rE^2*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*z*((3*J2*rE^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (mu*((18*J2*rE^2*y*z^3)/(x^2 + y^2 + z^2)^4 - (J2*rE^2*y*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^2 - (6*J2*rE^2*y*z)/(x^2 + y^2 + z^2)^3 + (4*J2*rE^2*y*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (3*mu*y*z*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0;...
                (mu*((J2*rE^2*((12*x*z)/(x^2 + y^2 + z^2)^2 - (24*x*z^3)/(x^2 + y^2 + z^2)^3))/(2*(x^2 + y^2 + z^2)) - (6*J2*rE^2*x*z^3)/(x^2 + y^2 + z^2)^4 + (J2*rE^2*x*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^2 - (4*J2*rE^2*x*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*z*((3*J2*rE^2*x*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*x*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) + (mu*x*((J2*rE^2*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (J2*rE^2*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x*z*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), (mu*((J2*rE^2*((12*y*z)/(x^2 + y^2 + z^2)^2 - (24*y*z^3)/(x^2 + y^2 + z^2)^3))/(2*(x^2 + y^2 + z^2)) - (6*J2*rE^2*y*z^3)/(x^2 + y^2 + z^2)^4 + (J2*rE^2*y*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^2 - (4*J2*rE^2*y*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*z*((3*J2*rE^2*y*z^2)/(x^2 + y^2 + z^2)^3 + (J2*rE^2*y*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) + (mu*y*((J2*rE^2*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (J2*rE^2*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*y*z*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), (mu*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(3/2) + (mu*((J2*rE^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2 - (J2*rE^2*(6/(x^2 + y^2 + z^2) - (30*z^2)/(x^2 + y^2 + z^2)^2 + (24*z^4)/(x^2 + y^2 + z^2)^3))/(2*(x^2 + y^2 + z^2)) + (2*J2*rE^2*z*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^2 - (4*J2*rE^2*z^2*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) + (2*mu*z*((J2*rE^2*((6*z)/(x^2 + y^2 + z^2) - (6*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)) - (J2*rE^2*z*((3*z^2)/(x^2 + y^2 + z^2) - 1))/(x^2 + y^2 + z^2)^2))/(x^2 + y^2 + z^2)^(3/2) - (3*mu*z^2*((J2*((3*z^2)/(x^2 + y^2 + z^2) - 1)*rE^2)/(2*x^2 + 2*y^2 + 2*z^2) - 1))/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0];
        end
        
        function H = getH(obj, state, lla, time)
            [H, ~] = gnc.ekf_update_mat(state, lla, time );
            
        end
    
    
    end
    

    
    
    
    
end