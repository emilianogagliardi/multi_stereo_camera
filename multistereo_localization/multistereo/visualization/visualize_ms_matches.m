function fs = visualize_ms_matches(ims1, ims2, k1, k2, idx1, idx2, show_style)

if nargin < 7
    show_style = 'montage';
end

if size(k1, 1) == 3
    k1 = from_homogeneous(k1);
end
if size(k2, 1) == 3
    k2 = from_homogeneous(k2);
end
fs = [];
for ii = 1:size(ims1, 3)
	for jj = 1:size(ims2, 3)
        idxii = idx1 == ii;
        idxjj = idx2 == jj;
        idx = idxii & idxjj;
        if any(idx) > 0
            matched1 = k1(:, idx);
            matched2 = k2(:, idx);
            fs = [fs, figure()];
            showMatchedFeatures(ims1(:, :, ii), ims2(:, :, jj), matched1', matched2', show_style);
        end
    end
end
end

