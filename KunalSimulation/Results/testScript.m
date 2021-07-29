protons = [1 2 3; 1 2 4]'

allowed_x_start = [0.1 0.3 0.5 0.7 0.9];
allowed_x_end = allowed_x_start + 0.15;

allowed_y_start = allowed_x_start + 1;
allowed_y_end = allowed_y_start + 0.15;

protons(1,:);

allowed = zeros(1, size(protons,2));

for c = 1:size(protons, 2)
    allowed(c) = any(protons(1,c)>allowed_x_start & protons(1,c)<allowed_x_end) ...
                && any(protons(2,c)>allowed_y_start & protons(2,c)<allowed_y_end);
end

mask = allowed == 1

protons = protons(:, mask)


A = [0 1 2; 5 1 2; 6 1 2]

Y=prctile(A, 75)

sum(A>Y, 1)