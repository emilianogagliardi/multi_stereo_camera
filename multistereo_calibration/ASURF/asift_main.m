close all;
clear all;
clc;

% ALGORITHM OUTLINE

% TAKE ONLY THE FIRST HALF OF BOTH IMAGES X
% ELIMINATE THE POINTS ON THE BOUNDARY X
% ELIMINATE THE MATCHES THAT ARE CAUSED BY SIMPLE VERTICAL MIRRORING
% TAKE THE BEST MATCHES?
% EMPLOY SOME COLOR DECISION ALGORITHM?
% EMPLOY SOME HOUGH TRANSFORM ON INTERSECTION WITH LICENSE PLATE?


% Algorithms params
METHOD = "SURF"; % SURF
FRAC = 1;

CLEAN = true;
FP  = 0.1;

SAME_T = 150;

im1_original = imread('./samples/evoque-back.jpg');
im1 = im1_original;

% pre) create im2_original flipping im1_original

[im2_original,A,Ai] = flipImage(im1_original);
im2 = im2_original;

[ro,co,~] = size(im1);

%1) taking only a relevant fraction of both images

im1 = im1(:,1:round(FRAC*co),:);
im2 = im2(:,1:round(FRAC*co),:);

r = ro;
c = round(FRAC*co);


tic
fprintf('Computing ASURF features for the original image. \n');
[kps1,descrs1] = affineDetect(im1,true);
fprintf('For the original image %d points have been found.\n',length(kps1));
toc

tic
fprintf('Computing ASURF features for the mirrored image. \n');
[kps2,descrs2] = affineDetect(im2,true);
fprintf('For the mirrored half %d points have been found.\n',length(kps2));
toc

% CLEANING PHASE: removing all the points that are in the outer frame of
% the image (which percentage is described by FP)

if(CLEAN == true)
    [kps1_c,descrs1_c] = cleanPoints(kps1,descrs1,r,c,FP);
    fprintf('For the original half %f%% points have been removed, %d left.\n',...
        100*((length(kps1)-length(kps1_c)) / length(kps1)), length(kps1_c));

    [kps2_c,descrs2_c] = cleanPoints(kps2,descrs2,r,c,FP);
    fprintf('For the mirrored half %f%% points have been removed, %d left.\n',...
        100*((length(kps2)-length(kps2_c)) / length(kps2)), length(kps2_c));
else
    fprintf('No cleaning phase performed.\n');
    kps1_c = kps1;
    descrs1_c = descrs1;
    kps2_c = kps2;
    descrs2_c = descrs2;
end

% matching features brute force
if(strcmp(METHOD,'SURF') || strcmp(METHOD,'OSURF'))
    fprintf('Matching features. \n');
    [indexPairs,scores] = matchFeatures(descrs1_c,descrs2_c,'Unique',true,...
        'Method','Approximate','MatchThreshold',1.0);
elseif(strcmp(METHOD,'SIFT') || strcmp(METHOD,'PHOW'))
    % da and db are the descriptors as they come out from [fa, da] = vl_sift(Ia)
    [matches, scores] = vl_ubcmatch(descrs1_c',descrs2_c',2.5) ;
    indexPairs = matches';
end

p1 = kps1_c(indexPairs(:,1),:);
p2 = kps2_c(indexPairs(:,2),:);
fprintf('A total of %d correspondances have been found.\n',length(indexPairs));

% remove matches that when backpropagated leads to the same point
fprintf('Removing useless matches. \n');
p1_final = [];
p2_final = [];
scores_ok = [];

for i=1:length(p1)
    p1o = p1(i,:);
    p2m = p2(i,:);
    p2bp = Ai*[p2m(1:2),1]';
    if(get2DDistance(p1o,p2bp(1:2)) > SAME_T)
        p1_final(end+1,:) = p1o;
        p2_final(end+1,:) = p2m;
        scores_ok(end+1,:) = scores(i);
    end
end

% extract the RANSAC homography
%{
[H,inl] = ransacfithomography(p1', p2', 0.01);

% extracting good ones
p1_ok = p1(inl,:);
p2_ok = p2(inl,:);
%}

% extracting best scores
%{
[~,idx] = sort(scores_ok);
TOP = min(length(idx),50);
topIdx = idx(1:TOP);
p1_ok = p1_final(topIdx,:);
p2_ok = p2_final(topIdx,:);
%}

p1_ok = p1_final;
p2_ok = p2_final;

%{
% showing correspondances
I1 = rgb2gray(im1);
I2 = rgb2gray(im2);

I = zeros([r c*2 1]);
I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
figure, imshow(I/255); hold on;

for i=1:length(p1_ok)
  cr=rand(1,3);
  %if( p1_ok(i,1) >= 12 && p1_ok(i,1) <= 112 && ...
          %p1_ok(i,2) >= 112 && p1_ok(i,2) <= 206)
      plot([p1_ok(i,1) p2_ok(i,1)+size(I1,2)],[p1_ok(i,2) p2_ok(i,2)],'-','Color',cr)
      plot([p1_ok(i,1) p2_ok(i,1)+size(I1,2)],[p1_ok(i,2) p2_ok(i,2)],'o','Color',cr)
  %end
end
%}

% take back the point of the mirrored image to their original domain

As = p1_ok;
Bs = zeros(length(As),2);

for i=1:length(p2_ok)
    B = Ai*[p2_ok(i,1:2),1]';
    Bs(i,:) = B(1:2);
end

figure;
imshow(im1_original);
hold on;

%{
for i=1:length(As)
    cr=rand(1,3);
    plot([As(i,1) Bs(i,1)],[As(i,2) Bs(i,2)],'-','Color',cr)
    plot([As(i,1) Bs(i,1)],[As(i,2) Bs(i,2)],'o','Color',cr)
end
%}


fprintf('Detect actual simmetries using RANSAC+Hough Transform. \n');

[H,inl] = ransacfithomography(As', Bs', 0.1);

fprintf("RANSAC: %d inliers found.\n",length(inl));

lines = zeros(3,length(inl));
gAs = zeros(length(inl),2);
gBs = zeros(length(inl),2);

for idx=1:length(inl)
   
    i = inl(idx);
    
    A = [As(i,:),1];
    B = [Bs(i,:),1];
    gAs(idx,:) = A(1:2);
    gBs(idx,:) = B(1:2);
    
    l = cross(A,B);
    lines(:,idx) = l;
    
end
[v,inl] = findVP(lines);

% plotting the correct simmetries
for idx=1:length(inl)
    cr=rand(1,3);
    i = inl(idx);
    plot([gAs(i,1) gBs(i,1)],[gAs(i,2) gBs(i,2)],'-','Color',cr)
    plot([gAs(i,1) gBs(i,1)],[gAs(i,2) gBs(i,2)],'o','Color',cr)
end

fprintf("HOUGH: %d inliers found.\n",length(inl));
plot(v(1),v(2),'xr','LineWidth',3);
fprintf('Completed. \n');