%Function to return starting positions of clusters within a structure.

function clusters = InitialiseClusters_QuickDiffusion(num_of_mol, electrodeactivation_flag, castellation_flag, sensor_xsize, sensor_ysize, solution_height, castellation_trappingregion_index, chip, r_initial)

    clusters.centre_x = zeros(1, num_of_mol);
    clusters.centre_y = zeros(1, num_of_mol);
    clusters.centre_z = zeros(1, num_of_mol);
    clusters.radius = r_initial*ones(1, num_of_mol);
    clusters.radius_prev = zeros(1,num_of_mol);

    %If electrodeactivation_flag = 0 (electrodes off), place DNA randomly
    %within bulk of reaction chamber
    if electrodeactivation_flag == 0
        for i = 1:num_of_mol

            %Clusters can take any x, y, or z value
            clusters.centre_x(i) = sensor_xsize*rand;
            clusters.centre_y(i) = sensor_ysize*rand;
            clusters.centre_z(i) = solution_height*rand;
        end

    %If electrodeactivation_flag = 1 (electrodes on), place DNA randomly along
    %the edges of the line electrodes, to which they would be attracted
    elseif electrodeactivation_flag == 1 && castellation_flag == 1
        for i = 1:num_of_mol

            %Select random trapping region from array
            trapping_region = randi(length(castellation_trappingregion_index));

            %Set x and y coords of DNA cluster centre from castellation
            %trapping region
            clusters.centre_x(i) = castellation_trappingregion_index(1, trapping_region);
            clusters.centre_y(i) = castellation_trappingregion_index(2, trapping_region);
            
            %Clusters start on top of electrode
            clusters.centre_z(i) = chip.electrode_thickness;
            
        end
    else
        
        for i = 1:num_of_mol
            
            %Select a random electrode for cluster origin. Then randomly select
            %one of the two edges of the electrode for cluster origin.
            electrode_origin = randi(chip.N_electrodes) - 1;
            edge_origin = randi(2) - 1;

            %Compute x position. Chip layout goes in repeating units of
            %electrode_width, electrodeisfet_separation, isfet_width,
            %electrodeisfet_separation
            clusters.centre_x(i) = electrode_origin*(chip.electrode_width + chip.electrode_separation) + edge_origin*chip.electrode_width;

            %Clusters can take any y position
            clusters.centre_y(i) = sensor_ysize*rand;
            
            %Clusters start on top of electrode
            clusters.centre_z(i) = chip.electrode_thickness;
        end
    end
end