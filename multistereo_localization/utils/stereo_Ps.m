function [Pl, Pr] = stereo_Ps(stereo)
%STEREO_PS get stereo projection matrices in left camera frame
[Kl, Kr, Rs, ts] = convert_stereo_params(stereo);
Pl = [Kl, zeros(3, 1)];
Pr = Kr * [Rs, ts];
end

