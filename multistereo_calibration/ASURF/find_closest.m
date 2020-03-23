% which are the two closest points to the one i would select?
% which are the regular descriptor of one and the mirrored of the other?
% which is their distance?

xl = [208.9915, 250.4344];
xr = [531.3765, 223.2453];

% get the closest to xl
mind = inf;
minidx = -1;
cont = 0;
idxl = [];
for i=1:length(kp1)
    d = get2DDistance(xl,kp1(i,:));
    if(d<=mind)
        mind = d;
        minidx=i;
    end
    
    if(d<5)
        idxl(end+1) = i;
    end
    
end

xli = kp1(minidx,:);
dl = desc(minidx,:);

% get the closest to xr
mind = inf;
minidx = -1;
idxr = [];

for i=1:length(kp1)
    d = get2DDistance(xr,kp1(i,:));
    if(d<=mind)
        mind = d;
        minidx=i;
    end
    
    if(d<5)
        idxr(end+1) = i;
    end
end

xri = kp1(minidx,:);
dr = desc_m(minidx,:);

dists = [];
for i=1:length(idxl)
    for j=1:length(idxr)
        d1 = desc_m(idxl(i),:);
        d2 = desc(idxr(j),:);
        dists(end+1)=sum(sqrt((d1-d2).^2));
    end
end