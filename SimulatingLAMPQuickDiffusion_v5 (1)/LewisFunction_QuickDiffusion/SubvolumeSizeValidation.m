%Function to validate the input of subvolume_size against the size of the
%smallest feature of the chip, and ensures an integer number of subvolumes
%fit within the reaction chamber

function subvolume_size = SubvolumeSizeValidation (subvolume_size, chip, sensor_xsize, sensor_ysize, solution_height)

    chip_temp = struct2cell(chip);
    chip_min = min([chip_temp{:}]);
    subvolume_flag1 = 0;
    subvolume_flag2 = 0;

    while subvolume_flag1 ~=1 || subvolume_flag2 ~=1

        %Check that subvolume size is appropriately small for chip layout
        while subvolume_size < chip_min
            subvolume_size = input(sprintf('Subvolume size is larger than recommended (> smallest chip feature). Enter new value (< %d recommended) or re-enter old value:\n', param_min));
            subvolume_flag2 = 0;
        end
        subvolume_flag1 = 1;

        %Check that subvolume size produces integer number of subolumes
        while round(sensor_xsize/subvolume_size) ~= sensor_xsize/subvolume_size && round(sensor_ysize/subvolume_size) ~= sensor_ysize/subvolume_size && round(solution_height/subvolume_size) ~= solution_height/subvolume_size
            subvolume_size = input('Subvolume size does not fit an integer number of subvolumes within reaction space. Enter a new value:');
            subvolume_flag1 = 0;
        end
        subvolume_flag2 = 1;

    end
end
