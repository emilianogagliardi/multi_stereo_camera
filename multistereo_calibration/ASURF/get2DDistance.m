
function [d]=get2DDistance(a,b)
    d = sqrt((a(1)-b(1)).^2 + (a(2)-b(2)).^2);
end
