function color = color_map(i)
%COLOR_MAP transform integer to color
switch(i)
    case 1
        color = [1, 0, 0];
    case 2
        color = [0, 1, 0];
    case 3
        color = [0, 0, 1];
    case 4
        color = [0, 0, 0];
    otherwise
        color = rand(1, 3);
end
        
end

