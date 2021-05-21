tic;

clusters = clusters_init;

f3 = figure();

Z = zeros(N,N);

s = surf(Z);
colorbar;
xlim([1 N]);
xlabel('x', 'FontSize', 16)
ylim([1 N]);
ylabel('y', 'FontSize', 16)
zlim([-0.01 1.5*est_tot_protons/(N*N)]);
zlabel('# of protons', 'FontSize', 16)
set(gca,'FontSize',16);

mystr = {'Projection of Release of Protons onto 2D Plane',['t = 0 minutes']};
    title(mystr, 'FontSize', 20);

track = zeros(length(t), N, N);

persisted_protons = 0;

% Keep track of actual # of protons released
tot_protons = 0;

for i = 2:length(t)
    % For every cluster...
    num_of_active_sig = sum(i>sigmoid_start);
    
    % For every cluster...
    for j = 1:num_of_active_sig
        
        % If cluster is active then update propogation front
        prev_r_h = clusters(j).radius_h;
        prev_r_v = clusters(j).radius_v;
        
        clusters(j).radius_h = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, sensor_dim, N);
        clusters(j).radius_v = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, solution_height, M);
        
        % Compute position of randomly distributed protons around propogation front
        numtorelease = round(burst*(volume_of_sphere(clusters(j).radius_h)  - volume_of_sphere(prev_r_h)));
        
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
                            < clusters(k).radius_h^2, : ) = [];
            end
        end

        clusters(j).protons_x = additions(:,1)';
        clusters(j).protons_y = additions(:,2)';
        clusters(j).protons_z = additions(:,3)';
        persisted_protons = persisted_protons + length(additions(:,1));
    end 
    
    % Update sensor plot
    s.ZData = sumSensor2(clusters, N);

    drawnow limitrate;
    
    % Update count of real total # of protons
    tot_protons = sum(s.ZData(:));
    
    % Keep track of sensor array along time
    %track(i-sigmoid_start(1)+1,:,:) = s.ZData;
    track(i,:,:) = s.ZData;
    
    mystr = {'Projection of Release of Protons onto 2D Plane',['t = ' num2str(t(i)/60) ' minutes']};
    title(mystr, 'FontSize', 20);
    % mymov(i) = getframe(f3);
end

toc;
%
% v = VideoWriter('project2d.avi');
% open(v);
% writeVideo(v, mymov);
% close(v);