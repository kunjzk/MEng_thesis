
function [d_intersection, p_intersection] =  LinePlane_Intersection( plane_normalvector, plane_point, line_vector, line_point)
    d_intersection = dot((plane_point - line_point), plane_normalvector)/dot( line_vector, plane_normalvector);%distance along diffusion vector to intersection. = inf if no intersection
    p_intersection = line_point + d_intersection*line_vector; %point of intersection
end