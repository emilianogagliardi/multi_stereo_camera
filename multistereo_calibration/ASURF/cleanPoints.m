function [kpsc,descrsc] = cleanPoints(kps,descrs,row,col,fp)
    %kpsc = [];
    %descrsc = [];
    keep = zeros(length(kps),1);
    parfor i=1:length(kps)
        p = kps(i,:);
        if (insideFrame(p,row,col,fp) == false)
            keep(i) = 1;
            %kpsc = vertcat(kpsc,p);
            %descrsc = vertcat(descrsc,descrs(i,:));
        end
    end
    
    outside = find(keep);
    kpsc = kps(outside,:);
    descrsc = descrs(outside,:);
end