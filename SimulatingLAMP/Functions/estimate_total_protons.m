function tot_protons = estimate_total_protons(N, M, sensor_dim, solution_height, speed_of_front, t, burst)


tic;

sigmoid_start = [1];
num_of_mol = 1;
clusters(1).centre_x = 16;
clusters(1).centre_y = 16;
clusters(1).centre_z = 16;
clusters(1).radius_h = 0.01;
clusters(1).radius_v = BuMeters(BuMeters(0.01,0,sensor_dim,N), 1,solution_height,M);
clusters(1).protons_x = [];
clusters(1).protons_y = [];
clusters(1).protons_z = [];

persisted_protons = 0;

% Keep track of actual # of protons released
tot_protons = 0;

% Start from when the first proton is released
for i = 1:length(t)
    
    num_of_active_sig = sum(i>sigmoid_start);
    
    % For every cluster...
    for j = 1:num_of_active_sig
        % If cluster is active then update propogation front
        prev_r_h = clusters(j).radius_h;
        prev_r_v = clusters(j).radius_v;
            
        clusters(j).radius_h = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, sensor_dim, N);
        clusters(j).radius_v = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, solution_height, M);
        
        % CHANGE THIS TO INTEGRAL OVER LAST 'TIME_STEP' SECONDS!!
        % Compute position of randomly distributed protons around propogation front
        numtorelease = round(burst*(volume_of_sphere(clusters(j).radius_h)  - volume_of_sphere(prev_r_h)));
        
        if numtorelease > 0
            c = 2*rand(numtorelease,1)-1;
            lon=2*pi*rand(numtorelease,1);
            lat=acos(c);
            a=cos(lon).*sin(lat);
            b=sin(lon).*sin(lat);

            additions = [clusters(j).radius_h*a'+clusters(j).centre_x;clusters(j).radius_h*b'+clusters(j).centre_y;clusters(j).radius_v*c'+clusters(j).centre_z]';

            % Remove the computed proton positions outside the array
            additions(additions(:,1)>N+1,:) = [];
            additions(additions(:,1)<1,:) = [];
            additions(additions(:,2)>N+1,:) = [];
            additions(additions(:,2)<1,:) = [];
            additions(additions(:,3)>M+1,:) = [];
            additions(additions(:,3)<1,:) = [];

            for k = 1:num_of_mol
                if(k~=j)
                    additions( (additions(:,1) - clusters(k).centre_x).^2 ...
                        + (additions(:,2) - clusters(k).centre_y).^2 ...
                        + (additions(:,3) - clusters(k).centre_z).^2 ...
                        <= clusters(k).radius_h^2, : ) = [];
                end
            end

            %Persist the protons onto the array
            clusters(j).protons_x = additions(:,1)';
            clusters(j).protons_y = additions(:,2)';
            clusters(j).protons_z = additions(:,3)';

            persisted_protons = persisted_protons + length(additions(:,1));

        end

        % Update sensor plot
        Z = sumSensor2(clusters, N);

        % Update count of real total # of protons
        tot_protons = sum(Z(:)) + persisted_protons;
    end
end

toc;



end

