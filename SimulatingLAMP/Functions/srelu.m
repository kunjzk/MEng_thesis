% Shifted Rectified Linear Unit
% i.e. f(x) = max(0, m*(x-c))
%
function f_x = srelu(x, m, c)

f_x = max(0, m*(x-c));

end

