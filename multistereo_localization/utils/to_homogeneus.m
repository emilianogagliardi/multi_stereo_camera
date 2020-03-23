function ph = to_homogeneus(p)
%TOHOMOGENEUS
% assume points as columns
ph = [p; ones(1, size(p, 2))];
end

