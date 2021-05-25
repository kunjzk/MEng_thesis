function clusters = InitializeClusterHeight(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial, z_fraction)
    
    clusters.radius = r_initial*ones(1, num_of_mol);
    clusters.radius_prev = zeros(1,num_of_mol);
         
    %Create random coordinates to test agianst trapping criteria
    clusters.centre_x = sensor_xsize*ones(1,num_of_mol).*0.5;
    clusters.centre_y = sensor_ysize*ones(1,num_of_mol).*0.5;
    clusters.centre_z = solution_height*ones(1,num_of_mol).*z_fraction;

end
