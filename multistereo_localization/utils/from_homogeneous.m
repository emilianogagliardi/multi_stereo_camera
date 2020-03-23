function p = from_homogeneous(ph)
%FROMHOMOGENEOUS 
% assume points expressed as colums
p = ph(1:end-1, :) ./ ph(end, :);
end

