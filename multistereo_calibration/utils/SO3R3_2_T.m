function T = SO3R3_2_T(p)
    R = rodrigues(p(1:3));
    t = p(4:6);
    if size(t, 1) ~= 3
        t = t';
    end
    T = [R, t; zeros(1, 3), 1];
end

