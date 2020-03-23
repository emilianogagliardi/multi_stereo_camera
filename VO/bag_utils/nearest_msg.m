function msg = nearest_msg(bag, stamp, topic, search_interval_size)
%NEAREST_MSG Returns the nearest message in time to the given stamp
% Sometimes it is almost the nearest 
%   Parameters:
%   - bag: a ros bag
%   - stamp: a struct with Sec and Nsec fields
    t = rostime2sec(stamp);
    inc = search_interval_size/2;
    initial = t - inc;
    final = t + inc;
    sel = select(bag, 'Topic', topic, 'Time', [initial, final]);
    msgs = readMessages(sel);
    msg = nan;
    dist = Inf;
    for ii = 1:length(msgs)
        t_ = rostime2sec(msgs{ii}.Header.Stamp);
        d = abs(t_ - t);
        if d < dist
            dist = d;
            msg = msgs{ii};
        end
    end
end
