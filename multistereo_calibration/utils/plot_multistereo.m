function plot_multistereo(Ts, stereos)
    figure, hold on;
    for ii = 1:size(Ts, 3)
        disp('hey');
        Rl = Ts(1:3, 1:3, ii);
        tl = Ts(1:3, 4, ii);
        [~, ~, Rs, ts] = convert_stereo_params(stereos{ii});
        Tr = [Rs, ts; zeros(1, 3), 1] * Ts(:, :, ii);
        Rr = Tr(1:3, 1:3);
        tr = Tr(1:3, 4);
        plot_camera(Rl, tl, 50, ['cam', num2str(ii), 'L'], [1, 0, 0]);
        plot_camera(Rr, tr, 50, ['cam', num2str(ii), 'R'], [0, 0, 1]);
        %plotCamera('Orientation', Rl, 'Location', -Rl'*tl, 'Size', 10, 'AxesVisible', true, 'Label', [num2str(ii), 'L']);
        %plotCamera('Orientation', Rr, 'Location', -Rr'*tr, 'Size', 10, 'Color', [0, 0, 1], 'AxesVisible', true, 'Label', [num2str(ii), 'R']);
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    axis tight;
end
