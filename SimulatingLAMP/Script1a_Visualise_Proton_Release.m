tic;

clusters = clusters_init;

% Create figure
f = figure;

Z = zeros(N, N);
track = zeros(length(t), N, N);

persisted_protons = 0;

% Keep track of actual # of protons released
tot_protons = 0;

s = surf(zeros(N,N));
xlim([1 N]);
ylim([1 N]);
zlim([1 N]);
hold on;

set(gca,'FontSize',16);
mystr = {'Release of Protons in 3D Space',['t = 0 minutes']};
    title(mystr, 'FontSize', 20);

xlabel('x', 'FontSize', 16)
ylabel('y', 'FontSize', 16)
zlabel('z', 'FontSize', 16)

for i = 2:length(t)
    
    num_of_active_sig = sum(i>sigmoid_start);
    cla(gca);
    
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
                    <= clusters(k).radius_h^2, : ) = [];
            end
        end
        
        
        y = scatter3(additions(:,1), additions(:,2), additions(:,3),'.', 'linewidth',2);
        
        %Persist the protons onto the array
        clusters(j).protons_x = additions(:,1)';
        clusters(j).protons_y = additions(:,2)';
        clusters(j).protons_z = additions(:,3)';
        
        persisted_protons = persisted_protons + length(additions(:,1));
        
    end
    
    [X,Y] = meshgrid(1:N);
    h = surf(X,Y,ones(N), 'FaceColor',[204, 255, 230]/255);
    %set(h,'facealpha',0.8);
    box on;
    
    drawnow limitrate;
    
    % Update sensor plot
    Z = sumSensor2(clusters, N);
    
    % Update count of real total # of protons
    tot_protons = sum(Z(:)) + persisted_protons;
    
    % Keep track of sensor array along time
    track(i,:,:) = Z;
    
    100*i/length(t);
    
    mystr = {'Release of Protons in 3D Space',['t = ' num2str(t(i)/60) ' minutes']};
    title(mystr, 'FontSize', 20);
    
    %mymov(i) = getframe(f);
end

toc;
% 
% %%
% v = VideoWriter('multicluster.avi');
% open(v);
% writeVideo(v, mymov);
% close(v);