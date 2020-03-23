function [Iflip,A,Ai] = flipImage(I)
%FLIPIMAGE Flips the image passed as an argument around the axis parallel
%to y and center in the centroid.

refAngle = 0;
p = [size(I,2)/2; size(I,1)/2]; %center of mass of the image
N = [cos(refAngle); sin(refAngle)]; % is this vector v?
d = dot(p,N);
A = [eye(2)-2*(N*N') 2*d*N; 0 0 1];
Ai = inv(A);

tform = affine2d(A');
xWorldLimits = [1 size(I,2)];
yWorldLimits = [1 size(I,1)];
J = imwarp(I,tform,'OutputView',imref2d(size(I),xWorldLimits,yWorldLimits));

Iflip = J;
end

