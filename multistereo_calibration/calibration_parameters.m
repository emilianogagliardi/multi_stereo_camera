
% A single image taken by a stereo camera should be the horizontal
% concatenation of the left and right image. Images must be named as 
% "camid-timestamp.png"
% where both camid and timestamp are numbers. If cam1 takes a picture
% of the pattern at time 10: 1-10.png
% The position of all the left camera will be expressed with respect to the
% left cam0
global IMAGES_FOLDER;
IMAGES_FOLDER = '/home/emiliano/Desktop/guidance/superdataset';
%IMAGES_FOLDER = '~/Development/matlab/multi_stereo_camera/multistereo_calibration/test/data_test_calibgraph';

global TEMPLATE;
TEMPLATE = imread('pattern_surf.png');
TEMPLATE = imresize(TEMPLATE, 0.05);

global TEMPLATE_SHORT_SIDE
TEMPLATE_SHORT_SIDE = 815.34; % mm


global MIN_INLIERS_SINGLE_IMAGE;
MIN_INLIERS_SINGLE_IMAGE = 25;

global MAX_REPR_ERR_GUESS;
MAX_REPR_ERR_GUESS = 1;
            
% keypoints extraction and matching
global AFFINE_INVARIANT_FEATURES;
AFFINE_INVARIANT_FEATURES = true;
global SURF_METRIC_THRESHOLD;
global NUM_OCTAVES_PATTERN;
global NUM_SCALE_LEVELS_PATTERN;
global MATCH_THRESHOLD;
global MAX_RATIO;
SURF_METRIC_THRESHOLD = 800;
NUM_OCTAVES_PATTERN = 3;
NUM_SCALE_LEVELS_PATTERN = 3; 
MATCH_THRESHOLD = 1;
MAX_RATIO = 0.6;

% % homography outliers rejection
% global HOMOG_MAX_DISTANCE;
% global HOMOG_TRIALS;
% HOMOG_MAX_DISTANCE = 1; % pixels
% HOMOG_TRIALS = 10000;

% p3p outliers rejection
global P3P_RANSAC_MAX_ITERATIONS;
global P3P_RANSAC_CONFIDENCE;
global P3P_RANSAC_MAX_REPROJ_ERROR;
P3P_RANSAC_MAX_ITERATIONS = 1000;
P3P_RANSAC_CONFIDENCE = 99;
P3P_RANSAC_MAX_REPROJ_ERROR = 1;

global OPTIMIZE_SCALE;
OPTIMIZE_SCALE = false;

global DEBUG_STEREO_GUESS;
DEBUG_STEREO_GUESS = true;