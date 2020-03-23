function p = T_2_SO3R3(T)
%T_2_SO3R3 Summary of this function goes here
p = [rodrigues(T(1:3, 1:3))', T(1:3, 4)'];
end

