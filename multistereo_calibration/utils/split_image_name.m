function [cam_idx, timestamp] = split_image_name(image_name)
    s = strsplit(image_name, '-');
    cam_idx = str2num(s{1});
    s = strsplit(s{2}, '.');
    timestamp = str2num(s{1});
end

