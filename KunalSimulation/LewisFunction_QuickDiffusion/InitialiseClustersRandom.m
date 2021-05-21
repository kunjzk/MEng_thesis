%Function to return starting positions of clusters within a structure.
%v2 allows setting a trapping reagion for DNA for more realistic trapping

function clusters = InitialiseClustersRandom(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial)

    clusters.radius = r_initial*ones(1, num_of_mol);
    clusters.radius_prev = zeros(1,num_of_mol);
         
    %Create random coordinates to test agianst trapping criteria
    clusters.centre_x = sensor_xsize*rand(1,num_of_mol);
    clusters.centre_y = sensor_ysize*rand(1,num_of_mol);
    clusters.centre_z = solution_height*rand(1,num_of_mol);
end