%% SimulatingLAMPKunal
% Moniri + Lewis with changes
% Castellations not considered
% For simplicity protons can pass through DNA clusters
% No Electrodes
% Trapping regions are used as before

%% Clear workspace and load functions
clear all; 
close all;

addpath('/Users/kunal/Desktop/thesis/KunalSimulation/Functions/')

time_steps = [1 3 5 10];

%% Flags

%Flag to plot initial distribution of clusters
plotflag_initial = 1;

%Flag for, and frequency of, plotting release of protons
plotflag_proton_release = 0;
plotflag_proton_frequency = 10;

plotflag_2d_surface = 0;
plot_2d_surface_frequency = 5;
save_2d_surface = 0;

stop_rxn_tolerance = 10;

%% Reaction
disp('.')

num_of_mol = 100; % # of initial molecules
speed_of_front = 1.73E-6 ; % metres per second
r_initial = 1.25E-6; % in metres

% Starting time to each cluster is Normal distribution. Values taken from
% fit to real data
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

removed_1 = zeros(1, 2401);
removed_3 = zeros(1, 801);
removed_5 = zeros(1, 481);
removed_10 = zeros(1,241);

not_removed_1 = zeros(1, 2401);
not_removed_3 = zeros(1, 801);
not_removed_5 = zeros(1, 481);
not_removed_10 = zeros(1,241);

released_1 = zeros(1, 2401);
released_3 = zeros(1, 801);
released_5 = zeros(1, 481);
released_10 = zeros(1, 241);

radius_1 = zeros(1, 2401);
radius_3 = zeros(1, 801);
radius_5 = zeros(1, 481);
radius_10 = zeros(1, 241);


saturation_1 = 0
saturation_3 = 0
saturation_5 = 0
saturation_10 = 0


for counter=1:4
    
    total_time = 40; % time in minutes
    time_step = time_steps(counter); % in seconds
    t = 0:time_step:total_time*60; % in seconds

    % Number of protons released is proportional to area of propagation front.
    % 'burst' is the proportionality constant. Can use estimateburst script to
    % estimate the parameter based on experimental observations.
    burst = 3E13; %1.15E9 for 6mV, 4E8 for pH

    % Create starting time dist (and sample from it)
    pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
    sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index 


    %% Diffusion

    %Effective 'rate constant' for H ion diffusion, k_diffusion.
    %Value from http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
    %for 50 degrees C. Units m^2/s.
    diffusivity_H = 1E-8;

    %% Chip design
    disp('.')

    %Design parameters
    %Line electrodes have width in x and length in y. Length is equal to the
    %sensor_ysize. For no electrodes, best to put electrode_width as zero and halve the
    %electrodeisfet_separation
    chip.N_x = 56;
    chip.N_y = 78; %N_x x N_y array

    % isfet array
    % Lets assume isfet separations are 1/5 of isfet dimension
    chip.isfetisfet_separation = 2E-6;
    % Distance between isfets and sensor (both dimensions)
    chip.isfet_startSeparation = 2E-6;
    chip.isfet_endSeparation = 2E-6;
    % ISFET dimensions
    chip.isfet_width = 35E-6;
    chip.isfet_length = 37E-6;
    % For plotting
    chip.height_unit = 20E-6;

    %Calculate the total x and y size of the sensor array
    sensor_xsize = chip.N_x*chip.isfet_width + (chip.N_x-1)*(chip.isfetisfet_separation) + chip.isfet_startSeparation + chip.isfet_endSeparation;
    sensor_ysize = chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.isfet_startSeparation + chip.isfet_endSeparation;

    %%

    %Separation between the edge of the sensor array and the wall of the reaction
    %chamber
    chip.wall_separation_xpos = 2E-6;
    chip.wall_separation_xneg = 2E-6;
    chip.wall_separation_ypos = 2E-6;
    chip.wall_separation_yneg = 2E-6;

    reaction_x_size = sensor_xsize + chip.wall_separation_xpos + chip.wall_separation_xneg;
    reaction_y_size = sensor_ysize + chip.wall_separation_ypos + chip.wall_separation_yneg;

    solution_vol = 12*10^(-9); %12uL
    solution_height = solution_vol/(sensor_xsize*sensor_ysize);
    sensor_vol = sensor_xsize*sensor_ysize*solution_height;

    % 1m^3 = 1000L
    avogadro = 6.02E23;
    initial_pH = 8.75;
    initial_conc = 10^(-1*initial_pH); %moles/L
    initial_n_prot = 10^(-1*initial_pH)*avogadro*solution_vol/1000; %molecules
    concentration = ((initial_n_prot/solution_vol)/avogadro)*1000; %moles/L
    initial_pH_reverse = log10(concentration)*-1;

    %% Detection Distance = Debye Length

    debye.perm_space = 8.8544E-12;
    debye.perm_water = 65;
    debye.k = 1.3806E-23;
    debye.T = 336.15;
    debye.N_a = 6.023E23;
    debye.e = 1.6022E-19;
    debye.I = 2;

    % Set the distance above ISFET sensors that capture/detection events will
    % occur
    debye_length = sqrt(debye.perm_space*debye.perm_water*debye.k*debye.T/(2*debye.N_a*debye.e^2*debye.I));
    detection_distance = debye_length;

    %% Permissable proton coordinates

    % Build set of allowable coordinates for protons to be sensed
    chip.x_allowed_start = zeros(1, chip.N_x);
    chip.x_allowed_end = zeros(1, chip.N_x);
    chip.x_allowed_start(1) = chip.wall_separation_xneg + chip.isfet_startSeparation;
    for c = 2:chip.N_x
        chip.x_allowed_start(c) = round(chip.x_allowed_start(c-1) + chip.isfet_width + chip.isfetisfet_separation, 10);
    end
    chip.x_allowed_end = chip.x_allowed_start + chip.isfet_width;

    chip.y_allowed_start = zeros(1, chip.N_y);
    chip.y_allowed_end = zeros(1, chip.N_y);
    chip.y_allowed_start(1) = chip.wall_separation_yneg + chip.isfet_startSeparation;
    for c = 2:chip.N_y
        chip.y_allowed_start(c) = chip.y_allowed_start(c-1) + chip.isfet_length + chip.isfetisfet_separation;
    end
    chip.y_allowed_end = chip.y_allowed_start + chip.isfet_length;

    %% End of Inputs

    %% Build Chip

    max_cluster_radius = sqrt(reaction_x_size^2 + reaction_y_size^2 + solution_height^2);

    %% Initialise Clusters

    %If DNA exists in solution, initialise starting positions of DNA clusters
    if num_of_mol ~= 0
        clusters = InitialiseClustersRandom(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial);
        %clusters = InitializeClusterHeight(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial, 0.5);
    end

    %Set start time of simulation. signmoid_start(1) to start at the first
    %amplification. Enter 1 to start at t=0.
    t_start = sigmoid_start(1);

    %% Initialise Sensors
    %Array to store the number of protons at each sensor over time. col =
    %sensor, row = time
    sensor_nH = zeros(chip.N_x*chip.N_y, length(t));
    cum_nH = zeros(chip.N_x*chip.N_y, length(t));
    total_protons_in_solution = zeros(1, length(t));

    %% Simulation
    % Array to store the x,y,z coords of all protons
    stop_rxn = 0;
    tic;
    steps_with_no_protons_released = 0;

    for i = t_start:length(t)-1
        
        total_prot_removed_in_this_step = 0;
        total_prot_free = 0;
        total_prot_released_in_this_step = 0;

        if steps_with_no_protons_released == stop_rxn_tolerance
            stop_rxn = 1;
        end

        if stop_rxn == 1
            disp('stop rxn flag')
            break
        end

        protons = [];
        %Carry out next step of LAMP reaction
        num_of_active_sig = sum(i>sigmoid_start);

        for j = 1:num_of_active_sig

            % If the radius of any cluster exceeds the maximum permissiable
            % radius (r > len(max dimension)/2), then the cluster is larger
            % than the solution. All reactants have been used up and protons
            % can no longer be produced. This shaves off alot of computation.
            if sum(clusters.radius > max_cluster_radius) > 0
                stop_rxn = 1;
            end

            if stop_rxn == 1
                disp('stop rxn flag')
                break
            end

            % If cluster is active then update propogation front
            clusters.radius_prev(j) = clusters.radius(j);
            clusters.radius(j) = speed_of_front*(t(i)-t(sigmoid_start(j))); 

            % Compute position of randomly distributed protons around propogation front
            numtorelease = round(burst*(volume_of_sphere(clusters.radius(j))  - volume_of_sphere(clusters.radius_prev(j))));
            %numtorelease = round(burst*(clusters.radius(j)^2));

            c = 2*rand(numtorelease,1)-1;
            lon=2*pi*rand(numtorelease,1);
            lat=acos(c);
            a=cos(lon).*sin(lat); %random points on a sphere
            b=sin(lon).*sin(lat);

            % These are the produced protons
            additions = [clusters.radius(j)*a'+clusters.centre_x(j);clusters.radius(j)*b'+clusters.centre_y(j); clusters.radius(j)*c'+clusters.centre_z(j)];
            total_released = size(additions,2);
            % Remove the computed proton positions outside the sensor
            additions(:,additions(1,:)>sensor_xsize+chip.wall_separation_xneg) = []; %this should be xneg
            additions(:,additions(1,:)<chip.wall_separation_xneg) = [];
            additions(:,additions(2,:)>sensor_ysize+chip.wall_separation_yneg) = [];
            additions(:,additions(2,:)<chip.wall_separation_yneg) = [];
            additions(:,additions(3,:)>solution_height) = [];
            additions(:,additions(3,:)<0) = [];
            

            % Remove protons that fall within other clusters
            for k = 1:num_of_mol
                if(k~=j)
                    additions(:,(additions(1,:) - clusters.centre_x(k)).^2 ...
                        + (additions(2,:) - clusters.centre_y(k)).^2 ...
                        + (additions(3,:) - clusters.centre_z(k)).^2 ...
                        <= clusters.radius(k)^2) = [];
                end
            end
            
            total_remaining = size(additions,2);

            protons = [protons additions];
            total_prot_removed_in_this_step = total_prot_removed_in_this_step + total_released-total_remaining;
            total_prot_free = total_prot_free + total_remaining;
            total_prot_released_in_this_step = total_prot_released_in_this_step + total_released;
        end

        if size(protons, 2) > 0

            steps_with_no_protons_released = 0;

            disp('protons released, i=')
            disp(i)
            % Diffusion Equation: Z(t+del_t) = Z(t) + sqrt(2.D.del_t).N(0,1)
            % Rearrange to get del_t, time taken to reach detection height
            del_t = (((detection_distance*ones(1, size(protons,2)) - protons(3,:))./normrnd(0,1, [1, size(protons,2)])).^2)./(2*diffusivity_H);

            % Use X, Y diffusion equations to compute final locations of protons
            % protons_final has 4 rows: x,y,z,t
            protons_final = zeros(size(protons, 1), size(protons, 2));
            protons_final(1,:) = protons(1,:) + sqrt(del_t).*(sqrt(2*diffusivity_H)*ones(1, size(protons,2))).*normrnd(0,1, [1, size(protons,2)]);
            protons_final(2,:) = protons(2,:) + sqrt(del_t).*(sqrt(2*diffusivity_H)*ones(1, size(protons,2))).*normrnd(0,1, [1, size(protons,2)]);
            protons_final(3,:) = detection_distance*ones(1, size(protons,2));
            protons_final(4,:) = del_t;
    %         disp('protons_final before removal:')
    %         disp(protons_final)

            % Remove protons outside permissable set of coordinates
            % If a proton hits a wall, remove it. Otherwise we need to consider
            % rebound angles

            protons_final = checkCoordinatesAndTime(protons_final, chip, time_step, i, length(t));

            % Figure out which sensor the protons landed on
            % Needs to be 1-indexed, use ceil instead of floow
            protons_final(1,:) = ceil((protons_final(1,:)-(chip.isfet_startSeparation+chip.wall_separation_xneg))./(chip.isfet_width+chip.isfetisfet_separation));
            protons_final(2,:) = ceil((protons_final(2,:)-(chip.isfet_startSeparation+chip.wall_separation_yneg))./(chip.isfet_length+chip.isfetisfet_separation));
            protons_final(4,:) = ceil(protons_final(4,:)./time_step);

            total_protons_in_solution(i) = total_protons_in_solution(i-1) + size(protons, 2)-size(protons_final,2);


            % Add protons to that sensor at the right time step     
    %             x = protons_final(1,c)
    %             y = protons_final(2,c)
    %             index = x+(y-1)*chip.N_x
    %             timesteps = i + protons_final(4,c):length(t) (protons persist
    %             from the time they were released till the end of the
    %             reaction)
            for c = 1:size(protons_final, 2)
                % Instantaneous
                sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c)) = ...
                    sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c)) + 1;

                % Cumulative
                cum_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c):length(t)) = ...
                    cum_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c):length(t)) + 1;
            end

        else
            disp('zero protons, i=')
            disp(i)
            total_protons_in_solution(i) = total_protons_in_solution(i-1);
            %If (1) there is a positive number of protons in the solution 
            % and (2) all clusters have started amplifying,
            % AND 
            % No new protons were released, this is a signal to stop rxn
            if max(total_protons_in_solution) > 0 && min(clusters.radius)>r_initial
               steps_with_no_protons_released = steps_with_no_protons_released + 1;
            end
        end
        
        %removed_prot_per_step = [removed_prot_per_step total_prot_removed_in_this_step];
        
        if counter == 1
            removed_1(i) = total_prot_removed_in_this_step;
            not_removed_1(i) = total_prot_free;
            released_1(i) = total_prot_released_in_this_step;
            %radius_1(i) = clusters.radius;
        end
        
        if counter == 2
            removed_3(i) = total_prot_removed_in_this_step;
            not_removed_3(i) = total_prot_free;
            released_3(i) = total_prot_released_in_this_step;
            %radius_3(i) = clusters.radius;
        end
        
        if counter == 3
            removed_5(i) = total_prot_removed_in_this_step;
            not_removed_5(i) = total_prot_free;
            released_5(i) = total_prot_released_in_this_step;
            %radius_5(i) = clusters.radius;
        end
        
        if counter == 4
            removed_10(i) = total_prot_removed_in_this_step;
            not_removed_10(i) = total_prot_free;
            released_10(i) = total_prot_released_in_this_step;
            %radius_10(i) = clusters.radius;
        end
        
    end

    disp('end of reaction')
    disp(toc)
    
    if counter == 1
        t1 = t;
        cum_nH_1 = cum_nH;
        saturation_1 = i;
        %removed_1 = removed_prot_per_step;
    end
    
    if counter == 2
        t3 = t;
        cum_nH_3 = cum_nH;
        saturation_3 = i;
        %removed_3= removed_prot_per_step;
    end
    
    if counter == 3
        t5 = t;
        cum_nH_5 = cum_nH;
        saturation_5 = i;
        %removed_5 = removed_prot_per_step;
    end
    
    if counter == 4
        t10 = t;
        cum_nH_10 = cum_nH;
        saturation_10 = i;
        %removed_10 = removed_prot_per_step;
    end
    
    if counter == 5
        t0point1 = t;
        cum_nH_0point1 = cum_nH;
    end
end

%% Plot temporal to figure out region of interest
figure
plot(t1, sum(cum_nH_1), 'DisplayName', 'Time step = 1', 'LineWidth', 2)
hold on
plot(t3, sum(cum_nH_3), 'DisplayName', 'Time step = 3', 'LineWidth', 2)
plot(t5, sum(cum_nH_5), 'DisplayName', 'Time step = 5', 'LineWidth', 2)
plot(t10, sum(cum_nH_10), 'DisplayName', 'Time step = 10', 'LineWidth', 2)
hold off
grid on 
grid minor
legend
title('Temporal Signal vs Time, Varying Time Step')
xlabel('Time (s)')
ylabel('Number of Protons Detected by ISFET Array')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(0:1:length(removed_1)-1, removed_1, 'DisplayName', 'Time step = 1', 'LineWidth', 1)
hold on
plot(0:1:length(removed_3)-1, removed_3, 'DisplayName', 'Time step = 3', 'LineWidth', 1)
plot(0:1:length(removed_5)-1, removed_5, 'DisplayName', 'Time step = 5', 'LineWidth', 1)
plot(0:1:length(removed_10)-1, removed_10, 'DisplayName', 'Time step = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Number of protons removed per time step')
xlabel('Step number')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(0:1:length(not_removed_1)-1, not_removed_1, 'DisplayName', 'Time step = 1', 'LineWidth', 1)
hold on
plot(0:1:length(not_removed_3)-1, (not_removed_3), 'DisplayName', 'Time step = 3', 'LineWidth', 1)
plot(0:1:length(not_removed_5)-1, (not_removed_5), 'DisplayName', 'Time step = 5', 'LineWidth', 1)
plot(0:1:length(not_removed_10)-1, (not_removed_10), 'DisplayName', 'Time step = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Number of protons not removed per time step (released-removed)')
xlabel('Step number')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(t1, cumsum(not_removed_1), 'DisplayName', 'Time step = 1', 'LineWidth', 1)
hold on
plot(t3, cumsum(not_removed_3), 'DisplayName', 'Time step = 3', 'LineWidth', 1)
plot(t5, cumsum(not_removed_5), 'DisplayName', 'Time step = 5', 'LineWidth', 1)
plot(t10, cumsum(not_removed_10), 'DisplayName', 'Time step = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Cumulative Number of Protons Remaining (released-removed)')
xlabel('Time (s)')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)


%% 

figure
plot(0:1:saturation_1, cumsum(released_1(1:saturation_1+1)), 'DisplayName', 'Time step = 1', 'LineWidth', 1)
hold on
plot(0:3:saturation_3*3, cumsum(released_3(1:saturation_3+1)), 'DisplayName', 'Time step = 3', 'LineWidth', 1)
plot(0:5:saturation_5*5, cumsum(released_5(1:saturation_5+1)), 'DisplayName', 'Time step = 5', 'LineWidth', 1)
plot(0:10:saturation_10*10, cumsum(released_10(1:saturation_10+1)), 'DisplayName', 'Time step = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Cumulative Number of Protons Released by Cluster')
xlabel('Time (s)')
ylabel('Number of Protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%% 

figure
plot(t1, radius_1, 'DisplayName', 'Time step = 1', 'LineWidth', 1)
hold on
plot(0:1:length(radius_3)-1, radius_3, 'DisplayName', 'Time step = 3', 'LineWidth', 1)
plot(0:1:length(radius_5)-1, radius_5, 'DisplayName', 'Time step = 5', 'LineWidth', 1)
plot(0:1:length(radius_10)-1, radius_10, 'DisplayName', 'Time step = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Cumulative Number of Protons Released by Clusters')
xlabel('Time (s)')
ylabel('Number of Protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

mean(removed_1)
mean(removed_3)
mean(removed_5)
mean(removed_10)
