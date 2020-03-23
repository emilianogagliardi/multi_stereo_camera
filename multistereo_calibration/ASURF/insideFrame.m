function [f]=insideFrame(p,row,col,fp)
    % point p expressed as (x,y)
    xbl = floor(col*fp);
    xbr = floor(col*(1-fp));
    ybt = floor(row*fp);
    ybb = floor(row*(1-fp));
    
    if(p(1) <= xbl || p(1) >= xbr || ...
            p(2) <= ybt || p(2) >= ybb)
        f = true;
        return;
    end
    
    f = false;
end