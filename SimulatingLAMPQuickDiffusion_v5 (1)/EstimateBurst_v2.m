%Estimate burst parameter
%Estimate the surface area of DNA clusters by counting the number of
%subvolumes surrounding each subvolume within a DNA cluster that is not
%within a DNA cluster, at each time step. From this, give an estimate for the
%proportionality constant that governs the number of H ions released due to
%a given increase in surface area of a DNA cluster as it grows during a
%LAMP reaction in order to meet experimentally observed pH decrease.
%v2 speeds up program by restructuring clusters.

clear

num_of_mol = 1000; % # of initial molecules
solution_height = 2.5E-3; % in meters
sensor_xsize = 4E-3;
sensor_ysize = 4E-3;
speed_of_front = 1.73E-6 ; % meters per second

clusters.centre_x = zeros(1, num_of_mol);
clusters.centre_y = zeros(1, num_of_mol);
clusters.centre_z = zeros(1, num_of_mol);
clusters.radius = 1.25E-6*ones(1, num_of_mol);

for i = 1:num_of_mol

    %Clusters can take any x, y, or z value
    clusters.centre_x(i) = sensor_xsize*rand;
    clusters.centre_y(i) = sensor_ysize*rand;
    clusters.centre_z(i) = solution_height*rand;
end

total_time = 40; % time in minutes
time_step = 5; % in seconds
t = 0:time_step:total_time*60; % in seconds

% Starting time to each cluster is Normal distribution...
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

% Create starting time dist (and sample from it)
pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index 

pH_0 = 8.75;
pH_delta = 0.9;

buffer_capacity = 1;

%Initial # H ions due to pH of solution
N_A = 6.02214086E23;
nH_delta = (10^(-(pH_0 - pH_delta))-10^(-pH_0))*N_A*sensor_xsize*sensor_ysize*solution_height;

surface_area = zeros(1,length(t));
surface_area_change = zeros(1,length(t)-1);
surface_area(1,1) = 0;

subvolume_size = 5E-6;

subvolume_xyzsize = [subvolume_size, subvolume_size, subvolume_size];
%Note x and y swapped to follow meshgrid convention
subvolume_xyz = {subvolume_xyzsize(2)/2:subvolume_xyzsize(2):sensor_xsize-subvolume_xyzsize(2)/2, subvolume_xyzsize(1)/2:subvolume_xyzsize(1):sensor_ysize-subvolume_xyzsize(1)/2, subvolume_xyzsize(3)/2:subvolume_xyzsize(3):solution_height-subvolume_xyzsize(3)/2};
[subvolume_mesh3x, subvolume_mesh3y, subvolume_mesh3z] = meshgrid(subvolume_xyz{1,1}, subvolume_xyz{1,2}, subvolume_xyz{1,3});

%=1 for subvolume in DNA cluster, =2 for subvolume in DNA cluster known to
%be entirely surrounded by DNA cluster, 0 otherwise.
subvolume_incluster = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

    for i = sigmoid_start(1):length(t)
        
        disp(t(i))
        num_of_active_sig = sum(i>sigmoid_start);
        
        % For every cluster...
        for j = 1:num_of_active_sig

            % If cluster is active then update propogation front
            prev_r = clusters.radius(j);
            clusters.radius(j) = speed_of_front*(t(i)-t(sigmoid_start(j))); %BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, sensor_dim, N);
           
            %Calculate distance from each subvolume to cluster surface
            subvolumetoclusterradius_mag = ((subvolume_mesh3x - clusters.centre_x(j)).^2 +(subvolume_mesh3y - clusters.centre_y(j)).^2 +(subvolume_mesh3z - clusters.centre_z(j)).^2).^0.5 - clusters.radius(j);        

            %Find index of subvolumetoclusterradius values that are negative.
            %Note meshgrid swaps x and y dimensions
            subvolumeincluster_linear_temp = find(subvolumetoclusterradius_mag < 0);
            
            for k = subvolumeincluster_linear_temp
                if subvolume_incluster(k) == 0
                    subvolume_incluster(k) = 1;
                end
            end
        end
        
        % Find all subvolumes within clusters, except those that have
        % already been determiend to be entirely surrounded by DNA cluster
        subvolumeincluster_linear = find(subvolume_incluster == 1);
        
        if ~isempty(subvolumeincluster_linear)
            
            %Count the number of faces surrounding DNA clusters
            [subvolumeincluster_indexy, subvolumeincluster_indexx, subvolumeincluster_indexz] = ind2sub(size(subvolumetoclusterradius_mag), subvolumeincluster_linear);
            
            for j =1:length(subvolumeincluster_linear)
                
                flags = [0,0,0,0,0,0];
                
                %If neighbouring subvolume in specified direction is not
                %within a DNA cluster, the flag is set to 1
                
                %+y
                if subvolumeincluster_indexy(j) ~= size(subvolume_incluster,1) &&   subvolume_incluster(subvolumeincluster_indexy(j)+1, subvolumeincluster_indexx(j), subvolumeincluster_indexz(j)) == 0
                    flags(1) = 1;
                end
                
                %-y
                if subvolumeincluster_indexy(j) ~= 1 && subvolume_incluster(subvolumeincluster_indexy(j)-1, subvolumeincluster_indexx(j), subvolumeincluster_indexz(j)) == 0
                    flags(2) = 1;
                end      
                
                %+x
                if subvolumeincluster_indexx(j) ~= size(subvolume_incluster,2) && subvolume_incluster(subvolumeincluster_indexy(j), subvolumeincluster_indexx(j)+1, subvolumeincluster_indexz(j)) == 0
                    flags(3) = 1;
                end
                
                %-x
                if subvolumeincluster_indexx(j) ~= 1 && subvolume_incluster(subvolumeincluster_indexy(j), subvolumeincluster_indexx(j)-1, subvolumeincluster_indexz(j)) == 0
                    flags(4) = 1;
                end
                
                %+z
                if subvolumeincluster_indexz(j) ~= size(subvolume_incluster,3) && subvolume_incluster(subvolumeincluster_indexy(j), subvolumeincluster_indexx(j), subvolumeincluster_indexz(j)+1) == 0
                    flags(5) = 1;
                end
                
                %-z
                if subvolumeincluster_indexz(j) ~= 1 && subvolume_incluster(subvolumeincluster_indexy(j), subvolumeincluster_indexx(j), subvolumeincluster_indexz(j)-1) == 0
                    flags(6) = 1;
                end
                
                %Add the total number of neighbouring subvolumes not within
                %a DNA cluster to the surface area. Conversion to actual
                %area happens later
                surface_area(1,i) = surface_area(1,i) + sum(flags);
                
                %If subvolume entirely surrounded by cluster, remove from
                %future computation
                if isempty(find(flags == 1,1))
                    subvolume_incluster(subvolumeincluster_indexy(j), subvolumeincluster_indexx(j), subvolumeincluster_indexz(j)) = 2;
                end
            end
            
            %{
            subvolumeincluster_surround_linear = find(subvolume_incluster == 2);
            [subvolumeincluster_surroundy,subvolumeincluster_surroundx,subvolumeincluster_surroundz] = ind2sub(size(subvolumetoclusterradius_mag), subvolumeincluster_surround_linear);
            
            
            figure(1)
            hold on
            scatter3(subvolumeincluster_indexx,subvolumeincluster_indexy,subvolumeincluster_indexz)
            scatter3(subvolumeincluster_surroundx,subvolumeincluster_surroundy,subvolumeincluster_surroundz)
            hold off
            shg
            %}
        end
    end
    
%Calculate total change in surface area and estimate burst parameter    
surface_area = surface_area*subvolume_size^2;
for i = 2:length(t)
    surface_area_change(i-1) = surface_area(i) - surface_area(i-1);
end
surface_area_change_total = sum(surface_area_change);

burst_fit = nH_delta/surface_area_change_total;
fprintf('Estimated burst = %d', burst_fit)
