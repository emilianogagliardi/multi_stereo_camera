
%------------------------- EXTRACTION PARAMETERS --------------------------
global SURF_THRESH;
global SURF_OCTAVES;
global SURF_SCALES;
% boolean, if true keypoints are uniformly selected after their extraction
% in left images
global SELECT_UNIFORM;
% number of keypoints to keep if SELECT_UNIFORM is true
global KEYPOINTS_NUMBER;
SURF_THRESH = 200;
SURF_OCTAVES = 5;
SURF_SCALES = 6;
SELECT_UNIFORM = true;
KEYPOINTS_NUMBER = 70;

%------------------------ STEREO MATCHING PARAMETERS ----------------------
% max descriptor distance
global STEREO_DESC_DISTANCE;
% used only if USE_EPIPOLAR_SEARCH = false, discard matching between d1 and
% d2 if the distance between d1 and its second best match is grather than 
% STEREO_MATCH_THRESHOLD * norm(d1-d2)
global STEREO_MAX_RATIO;
% boolean, true if fast stereo matching on rectified images
global USE_EPIPOLAR_SEARCH;
% precision of epipolar search
global EPIPOLAR_TICKNESS;
% maximum allowed disparity between a stereo pair of keypoints
global MAX_STEREO_DISPARITY;
% boolean, remove bad stereo matching after triangulation. Remove all
% points in triangulation that have reprojection error higher than the mean
% reprojection error among all the triangulated points
global REPROJ_ERROR_OUTLIERS;

STEREO_DESC_DISTANCE = 1.2;
STEREO_MAX_RATIO = 0.8;
USE_EPIPOLAR_SEARCH = false;
EPIPOLAR_TICKNESS = 2;
MAX_STEREO_DISPARITY = inf;
REPROJ_ERROR_OUTLIERS = true;

%------------------- KEYPOINTS TRACKING PARAMETERS ------------------------
% window size of search in pixels for keyopints tracking by constant motion
% model (square side)
global TRACK_SEARCH_SIZE;
TRACK_SEARCH_SIZE = 20;

%---------- ABSOLUTE ORIENTATION OUTLIERS REJECTION PARAMETERS ------------
global ABSOR_RANSAC_SAMPLE_NUMBER;
global ABSOR_RANSAC_MAX_ITERATIONS;
global ABSOR_RANSAC_CONFIDENCE;
global ABSOR_RANSAC_MAX_3D_DISTANCE;
global ABSOR_RANSAC_MAX_REPROJ_ERR; % when repr err is used in outliers rej

ABSOR_RANSAC_SAMPLE_NUMBER = 3;
ABSOR_RANSAC_MAX_ITERATIONS = 1000;
ABSOR_RANSAC_CONFIDENCE = 99;
ABSOR_RANSAC_MAX_3D_DISTANCE = 300^2; % mm to the square
ABSOR_RANSAC_MAX_REPROJ_ERR = 2^2; % pixels to the square

%------------- SINGLE CAMERA OUTLIERS REJECTION (P3P) PARAMETERS ----------
global P3P_RANSAC_MAX_ITERATIONS;
global P3P_RANSAC_CONFIDENCE;
global P3P_RANSAC_MAX_REPROJ_ERROR;

P3P_RANSAC_MAX_ITERATIONS = 1000;
P3P_RANSAC_CONFIDENCE = 99;
P3P_RANSAC_MAX_REPROJ_ERROR = 8^2; % pixels to the square

%-------------- REPROJECTION ERROR MINIMIZATION PARAMETERS ----------------
global REPROJ_MINIMIZ_METHOD; % either LM or LHM, multicamera always uses LM
global LHM_EPSILON; % error threshold to stop
global LHM_TOLERANCE; % stop if abs(old_error-new_err)/old_err < LHM_TOLERANCE
global LHM_MAX_ITERATIONS;
global LM_TOLX; % see optimoptions matlab fun, 'TolX' and 'TolFun'
global LM_TOLFUN;
global LM_MAX_ITERATIONS;

REPROJ_MINIMIZ_METHOD = 'LHM';
LHM_EPSILON = 1e-6;
LHM_TOLERANCE = 1e-6;
LHM_MAX_ITERATIONS = 10;
LM_TOLX = 1e-6;
LM_TOLFUN = 1e-6;
LM_MAX_ITERATIONS = 30;

%----------- MULTICAMERA 3D2D OUTLIERS REJECTION (GP3P) PARAMETERS --------
global GP3P_RANSAC_MAX_ITERATIONS;
global GP3P_RANSAC_CONFIDENCE;
global GP3P_RANSAC_MAX_REPROJ_ERROR;

GP3P_RANSAC_MAX_ITERATIONS = 1000;
GP3P_RANSAC_CONFIDENCE = 99;
GP3P_RANSAC_MAX_REPROJ_ERROR = 0.0005; % collinearity max error (see opengv)



