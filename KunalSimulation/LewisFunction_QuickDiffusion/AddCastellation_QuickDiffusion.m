%Function to add castellations to chip design and provide the index of
%trapping locations

function [castellation_trappingregion_index, chip_layout] = AddCastellation_QuickDiffusion(x, y, chip, chip_layout)

    %Number of castellations in row is N_y*ncastellation_per pixel.
    %N_electrode*2-2 rows. *2 for castellations on both sides of electrodes
    %-2 for no rows on left at beginning and on right at end
    ncastellation_perrow = chip.N_y*chip.ncastellation_perpixel;
    castellation_trappingregion_index = zeros(2,(chip.N_electrodes*2-2)*ncastellation_perrow*2);

    %Calculate bounding limits of castellations on each side of electrode
    for i = 1:chip.N_electrodes
        
        %Calculate bounding limits of line electrodes in x direction only
        electrode_xstart = (i-1)*(chip.electrode_width + 2* chip.electrodeisfet_separation + chip.isfet_width + 2*chip.castellation_extension);
        
        %For the two sets of castellation per electrode
        for j = 1:2

            % Skip castellations that would extend beyond reaction chamber
            if ~(j == 1 && i == 1 || j == 2 && i == chip.N_x+1)

                %Find x bounds of castellation
                castellation_xstart = electrode_xstart + (j-1)*chip.electrode_width + (j-2)*chip.castellation_extension;
                castellation_xend = electrode_xstart + (j-1)*chip.electrode_width + (j-1)*chip.castellation_extension;

                %Find indices where x axis is between these values
                castellation_xindex_linear = find(x > castellation_xstart & x < castellation_xend);

                %For the sets of castellations at each pixel
                for k = 1:chip.N_y

                    %For each castellation per pixel
                    for l = 1:chip.ncastellation_perpixel

                        %Find y bounds of castellation
                        castellation_ystart = chip.sensor_startdisplacement + 0.5*chip.isfet_length + (k-1)*(chip.isfet_length + chip.isfetisfet_separation)...
                            - chip.castellation_width_perpixel*0.5 + (l-1)*(chip.castellation_width + chip.castellation_separation);
                        castellation_yend = chip.sensor_startdisplacement + 0.5*chip.isfet_length + (k-1)*(chip.isfet_length + chip.isfetisfet_separation)...
                            - chip.castellation_width_perpixel*0.5 + (l-1)*(chip.castellation_width + chip.castellation_separation) + chip.castellation_width;

                        %Find indices of subvolumes between bounds
                        castellation_yindex_linear = find(y > castellation_ystart & y < castellation_yend);

                        %Add castellations to chip layout
                        for m = castellation_xindex_linear
                            for n = castellation_yindex_linear
                                chip_layout(m,n,2) = 1;
                                chip_layout(m,n,1) = 0;
                                chip_layout(m,n,3) = 0;
                            end
                        end

                        %Save the relevant corners of the castellations as
                        %trapping points for DNA
                        if j == 1
                            
                            %For castellations on '-x' side of electrode
                            castellation_trappingregion_index(:,(i-2)*ncastellation_perrow*2 + chip.ncastellation_perpixel*(k-1) + (l-1)*2:(i-2)*ncastellation_perrow*2 ...
                            + chip.ncastellation_perpixel*(k-1) + (l-1)*2 + 1) = [castellation_xstart, castellation_xstart; castellation_ystart, castellation_yend];
                        else
                            
                            %For castellations on 'x' side of electrode
                            castellation_trappingregion_index(:,(chip.N_electrodes*2-4+i)*ncastellation_perrow*2 + chip.ncastellation_perpixel*(k-1) + (l-1)*2:(i-2)*ncastellation_perrow*2 ...
                            + chip.ncastellation_perpixel*(k-1) + (l-1)*2 + 1) = [castellation_xend, castellation_xend; castellation_ystart, castellation_yend];
                        end      
                    end
                end
            end
        end
    end
end