function [kLm, kRm, dLm, dRm] = epipolar_search(descL, descR, kptsL, kptsR, norm_function)
%EPIPOLAR_SEARCH computes matches in rectified stereo images
% implementation taken from PROSLAM with small modification to handle
% subpixel accuracy, and max allowed disparity. sfot comparison are used to
% allow epipolar tickness greater than 1
% Parameters
% - descL: NxX where N is the left keypoints number and X is the descriptor
% dimensionality, left keypoints associated descriptors
% - descR: MxX where M is the right keypoints number and X is the descriptor
% dimensionality, right keypoints associated descriptors
% - kptsL: Nx2 where N is the left keypoints number 
% - kptsR: Mx2 where M is the right keypoints number
% - norm_function: @(x, y) -> d where x and y are descriptors and d a
% scalar representing their distance

% KEYPOINTS: [COL, ROW]

global EPIPOLAR_TICKNESS;
global STEREO_DESC_DISTANCE;
global MAX_STEREO_DISPARITY;

% sort keypoints according to rows, then remap descriptors
idxL_map = quicksort(kptsL, @(x, i, j) compare(x, i, j, EPIPOLAR_TICKNESS));
idxR_map = quicksort(kptsR, @(x, i, j) compare(x, i, j, EPIPOLAR_TICKNESS));
kptsL = kptsL(idxL_map, :);
kptsR = kptsR(idxR_map, :);
descL = descL(idxL_map, :);
descR = descR(idxR_map, :);

% all keypoints are available for matching
left_available = true(size(kptsL, 1), 1);
right_available = true(size(kptsR, 1), 1);

matches = zeros(min(size(kptsL, 1), size(kptsR, 1)), 2);
matches_idx = 1;
    
idxR = 1;

% loop over left kpts
for idxL = 1:size(kptsL, 1)    

    % if current left not yet matched
    if left_available(idxL)

        % if checked all right kpts skip this left point
        if idxR == size(kptsR, 1) + 1
            break
        end

        % right keypoints are on lower row, skip left
        while smaller(kptsL(idxL, 2), kptsR(idxR, 2), EPIPOLAR_TICKNESS) ...
                || ~left_available(idxL)
            idxL = idxL + 1;
            if idxL == size(kptsL, 1) + 1
                break;
            end
        end


        % if checked all the left keypoints exit
        if idxL == size(kptsL, 1) + 1
            break;
        end

       % right keypoints are on upper row, skip them
       while greather(kptsL(idxL, 2), kptsR(idxR, 2), EPIPOLAR_TICKNESS)
           idxR = idxR + 1;
            if idxR == size(kptsR, 1) + 1
                break;
            end
       end

       if idxR == size(kptsR, 1) + 1
            break
       end

       % search for matching
       idxR_search = idxR;
       best_distance = STEREO_DESC_DISTANCE;
       idxR_best = 0;
       % scan epipolar line
       while equal(kptsL(idxL, 2), kptsR(idxR_search, 2),  EPIPOLAR_TICKNESS)

           % disparity stop condition
           if kptsL(idxL, 1) > kptsR(idxR_search, 1) + MAX_STEREO_DISPARITY || ...
                   kptsL(idxL, 1) < kptsR(idxR_search, 1)
               break;
           end

           % compute distance and eventually update best match
           if right_available(idxR_search)
               distance = norm_function(descL(idxL, :), descR(idxR_search, :));
               if distance <= best_distance
                   best_distance = distance;
                   idxR_best = idxR_search;
               end
           end

           % increment search index and check if finished
           idxR_search = idxR_search + 1;
           if idxR_search == size(kptsR, 1) + 1
               break;
           end
       end

       % check if something was found
       if idxR_best ~= 0
           % store new match
           matches(matches_idx, 1) = idxL;
           matches(matches_idx, 2) = idxR_best;
           matches_idx = matches_idx + 1;
           % keypoints no more available
           left_available(idxL) = false;
           right_available(idxR_best) = false;
           % restart from the next right keypoint
           idxR = idxR_best + 1;
       end

    end % if left available 
end % loop over left

% matches filled up to matches_idx
matches = matches(1:matches_idx - 1, :);
kLm = kptsL(matches(:, 1), :);
kRm = kptsR(matches(:, 2), :);
dLm = descL(matches(:, 1), :);
dRm = descR(matches(:, 2), :);
end

%--------------------------------------------------------------------------
% helper functions for sorting

function o = compare(x, i, j, precision)
    if greather(x(i, 2), x(j, 2), precision)
        o = 1;
    elseif smaller(x(i, 2), x(j, 2), precision)
        o = -1;
    else
        if greather(x(i, 1), x(j, 1), precisionz)
            o = 1;
        elseif smaller(x(i, 1), x(j, 1), precision)
            o = -1;
        else
            o = 0;
        end
    end
end

function o = equal(a, b, precision)
    o = abs(a-b) < precision/2;
end

function o = greather(a, b, precision)
    o = a-b > precision/2;
end

function o = smaller(a, b, precision)
    o = greather(b, a, precision);
end


