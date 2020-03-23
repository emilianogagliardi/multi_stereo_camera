function [Kl, Kr, Rs, ts] = convert_stereo_params(s)
    Kl = s.CameraParameters1.IntrinsicMatrix';
    Kr = s.CameraParameters2.IntrinsicMatrix';
    Rs = s.RotationOfCamera2';
    ts = s.TranslationOfCamera2';
end

