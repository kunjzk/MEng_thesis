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
zlim([-0.01 2*est_tot_protons/(N*N)]);
zlabel('# of protons', 'FontSize', 16)
set(gca,'FontSize',16);

mystr = {'Estimate of Array Output',['t = 0 minutes']};

title(mystr, 'FontSize', 20);

track = zeros(length(t), N, N);

persisted_protons = 0;
persisted_protons2 = 0;
saturated = 0;

% Keep track of actual # of protons released
tot_protons = 0;

for i = 2:length(t)
    
    % For every cluster...
    num_of_active_sig = sum(i>sigmoid_start);
    
    if(~saturated)
        persisted = 0;
        % For every cluster...
        for j = 1:num_of_active_sig
            % If cluster is active then update propogation front
            prev_r_h = clusters(j).radius_h;
            prev_r_v = clusters(j).radius_v;
            
            real_radius = speed_of_front*(t(i)-t(sigmoid_start(j)));
            clusters(j).radius_h = BuMeters(real_radius, 1, sensor_dim, N);
            clusters(j).radius_v = BuMeters(real_radius, 1, solution_height, M);
            
            if(clusters(j).radius_h < sqrt(2)*N && clusters(j).radius_v < sqrt(2)*M)
                
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
                
                for k = 1:num_of_active_sig
                    if(k~=j)
                        additions( (additions(:,1) - clusters(k).centre_x).^2 ...
                            + (additions(:,2) - clusters(k).centre_y).^2 ...
                            + (additions(:,3) - clusters(k).centre_z).^2 ...
                            < clusters(k).radius_h^2, : ) = [];
                    end
                end
                
                %Persist the protons onto the array
                clusters(j).protons_x = additions(:,1)';
                clusters(j).protons_y = additions(:,2)';
                clusters(j).protons_z = additions(:,3)';
                
                persisted = persisted + length(additions(:,1));
            else
                clusters(j).protons_x = [];
                clusters(j).protons_y = [];
                clusters(j).protons_z = [];
            end
        end
        
        % Update sensor plot
        % s.ZData = sumSensor2(clusters, N) + sumSensor(persisted_protons, N);
        rndProt = sqrt(persisted_protons/(N*N))*randn(N,N)+persisted_protons/(N*N);
        rndProt(rndProt<0) = 0;
        s.ZData = sumSensor2(clusters, N) + rndProt;
    else
        s.ZData = sqrt(persisted_protons/(N*N))*randn(N,N)+persisted_protons/(N*N);
    end
    
    drawnow limitrate;
    
    % Update count of real total # of protons
    tot_protons = sum(s.ZData(:));
    
    persisted_protons = persisted_protons + persisted;
    
    % Update count of real total # of protons
    tot_protons = sum(Z(:));
    
    if(tot_protons>1000 && persisted==0)
        tot_prot_count = tot_prot_count + 1;
        if(tot_prot_count>3)
            saturated = 1;
        end
    else
        tot_prot_count = 0;
    end
    
    % Keep track of sensor array along time
    %track(i-sigmoid_start(1)+1,:,:) = s.ZData;
    track(i,:,:) = s.ZData;
    
    mystr = {'Estimate of Array Output',['t = ' num2str(t(i)/60) ' minutes']};
    
    title(mystr, 'FontSize', 20);
    %
    if(ismember(t(i), [1:19]*60))
        %drawnow;
        %savefig([pwd '/Figures4Grant/N' num2str(num_of_mol) '_Chip_' num2str(t(i)/60) 'min.fig']);
        %print([pwd '/Figures4Grant/N' num2str(num_of_mol) '_Chip_' num2str(t(i)/60) 'min.png'],'-dpng')
        %waitforbuttonpress;
    end
    % mymov(i) = getframe(f3);
end

toc;

%
% v = VideoWriter('final.avi');
% open(v);
% writeVideo(v, mymov);
% close(v);