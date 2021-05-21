%Function to 'build' chip design using input parameters. The location of
%sensors and the trapping regions are output, along with a RGB image of the
%chip for use in figures
%v2 outputs coordinates for isfet row and columns separately
function [chip_layout, column_location, row_location, sensor_plotlocation] = BuildChip_QuickDiffusionv2(x, y, chip, resolution)


    %In chip_layout. 3rd dimension for
    %RGB colour values. Isfet = red (chip_layout(:,:,1)=1). Electrodes = green
    %(chip_layout(:,:,2)=1)
    chip_layout = ones(length(x), length(y), 3);

    %column_location denotes each sensor, columns denote coordinates
    %that bound sensing area: x_start, x_end
    column_location = zeros(chip.N_x, 2); 
    
    %row_location denotes each sensor, columns denote coordinates
    %that bound sensing area: y_start, y_end
    row_location = zeros(chip.N_y, 2); 
    
    %Determine which subvolumes on z=0 plane are above ISFET sensing regions.
    %Row of sensor_location denotes each sensor, column denotes the index of
    %each subvolume on the sensor.
    sensor_plotlocation = zeros(chip.N_x*chip.N_y, chip.isfet_length*chip.isfet_width/resolution^2); 
    
    %Determine which subvolumes on z=0 plane are above ISFET sensing regions.
    for i = 1:chip.N_x
        for j = 1:chip.N_y

            %Calculate the bounding limits of the sensor in terms of x and y
            %co-ordinates
            sensor_xstart = i*(chip.electrode_width + chip.electrodeisfet_separation + chip.castellation_extension)...
                + (i - 1)*(chip.isfet_width + chip.electrodeisfet_separation + chip.castellation_extension);
            sensor_xend = i*(chip.electrode_width + chip.electrodeisfet_separation + chip.isfet_width + chip.castellation_extension)...
                + (i - 1)*(chip.electrodeisfet_separation + chip.castellation_extension);
            sensor_ystart = chip.sensor_startdisplacement + (j - 1)*(chip.isfet_length + chip.isfetisfet_separation);
            sensor_yend = chip.sensor_startdisplacement + j*chip.isfet_length + (j - 1)*(chip.isfetisfet_separation);

            %Save coordinates
            column_location(i, 1) = sensor_xstart;
            column_location(i, 2) = sensor_xend;
            row_location(j, 1) = sensor_ystart;
            row_location(j, 2) = sensor_yend;
            
            %Find elements that fall within x and y bounds
            sensor_xindex_linear = find(x >= sensor_xstart & x <= sensor_xend);
            sensor_yindex_linear = find(y >= sensor_ystart & y <= sensor_yend);
            
            %Convert x, y co-ordinates to linear indices for easy storage
            for k = 1:length(sensor_xindex_linear)
                for l = 1:length(sensor_yindex_linear)
                    
                    %Save linear indices
                    sensor_plotlocation((i-1)*chip.N_y + j,(l-1)*length(sensor_xindex_linear) + k) = sub2ind([length(x), length(y)], sensor_xindex_linear(k), sensor_yindex_linear(l));

                    %Add to chip and sensor layout
                    chip_layout(sensor_xindex_linear(k), sensor_yindex_linear(l), 1) = 1;
                    chip_layout(sensor_xindex_linear(k), sensor_yindex_linear(l), 2) = 0;
                    chip_layout(sensor_xindex_linear(k), sensor_yindex_linear(l), 3) = 0;

                end
            end
        end
    end

    %Determine location of electrodes to create chip layout image for scatter
    %plots. 
    for i = 1:chip.N_electrodes

        %Calculate bounding limits of electrodes in x direction only
        electrode_xstart = (i-1)*(chip.electrode_width + 2*chip.electrodeisfet_separation + chip.isfet_width + 2*chip.castellation_extension);
        electrode_xend = i*chip.electrode_width + (i-1)*(2*chip.electrodeisfet_separation + chip.isfet_width + 2*chip.castellation_extension);

        %Find indices where x axis is between these values
        electrode_xindex_linear = find(x > electrode_xstart & x < electrode_xend);

        %Add to chip layout
        for j = electrode_xindex_linear
            chip_layout(j,:,2) = 1;
            chip_layout(j,:,1) = 0;
            chip_layout(j,:,3) = 0;
        end
    end  
end
