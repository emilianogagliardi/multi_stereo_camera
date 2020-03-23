clear
close all
clc

load('stereos.mat');
cams = {cam1, cam2, cam3, cam4};

G = CalibrationGraph(cams);
[G, sol] = G.Calibrate();

%%
G2 = G.RemoveBadPoints(1);
[G2, sol2] = G2.Calibrate();

%%
[~, camera_idx, image_idx, all_errs] = G2.ComputeReprojectionErrors();

%% plot camera 1 reprojection errors
errs1 = all_errs(camera_idx == 1);
histogram(errs1, 'BinWidth', 0.05)
xlabel('Reprojection error [pixels]');
ylabel('Occurrences');
xticks([0:0.1:1.2]);
%% plot camera 2 reprojection errors
errs2 = all_errs(camera_idx == 2);
histogram(errs2, 'BinWidth', 0.05);
xlabel('Reprojection error [pixels]');
ylabel('Occurrences');
xticks([0:0.1:1.2]);
%% plot camera 3 reprojection errors
errs3 = all_errs(camera_idx == 3);
histogram(errs3, 'BinWidth', 0.05);
xlabel('Reprojection error [pixels]');
ylabel('Occurrences');
xticks([0:0.1:1.2]);
%% plot camera 4 reprojection errors
errs4 = all_errs(camera_idx == 4);
histogram(errs4, 'BinWidth', 0.05)
xlabel('Reprojection error [pixels]');
ylabel('Occurrences');
xticks([0:0.1:1.2]);
%% plot number of inliers for image
num = [];
for ii = 1:max(image_idx)
    n = sum(image_idx == ii);
    %if n > 19
        num = [num, n];
    %end
end
histogram(num, 'BinWidth', 10);

%% plot calibration graph
h = plot(G2.G, 'g', 'EdgeColor', [193,191,181]/255, 'NodeColor', [27,73,101]/255, 'NodeLabel', {});%, 'MarkerSize', 20, 'LineWidth', 4);%, 'Layout', 'circle');
highlight(h,{'Cam0', 'Cam1', 'Cam2', 'Cam3'}, 'NodeColor', [176,46,12]/255, 'Marker', 'd', 'MarkerSize', 20);
%%
xd = get(h, 'XData');
yd = get(h, 'YData');
text(xd(1:4)-0.1, yd(1:4)-0.1, {'Stereo1', 'Stereo2', 'Stereo3', 'Stereo4'}, 'FontSize',18, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
%text(xd(4:6)-0.1, yd(4:6)-0.1, {'PatternPose1', 'PatternPose2', 'PatternPose3'}, 'FontSize',18, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
