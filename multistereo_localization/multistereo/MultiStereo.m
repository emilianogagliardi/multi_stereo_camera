classdef MultiStereo
    
    properties (Access = private)
        Kls % calibration matrices left
        Krs % calibration matrices right
        % Tii = [Rs(:, :, ii), ts(:, ii)] brings points from cam1 left 
        % to camii left frame
        Rs 
        ts
        % Tsii = [stereosR(:, :, ii), stereost(:, ii)] brings points from
        % left cam to rigth cam of the stereo pair ii
        stereosR 
        stereost 
        % projection matrices in camera 1 left frame
        Pls
        Prs
        
        N % number of cameras
    end
    
    methods
        
        function obj = MultiStereo(stereos, multicalib)
            %MULTISTEREO Construct an instance of this class
            % Parameters:
            % - stereos: cell array of N stereoParameters (stereo camera
            % calibration matlab app). stereos(1) is considered as
            % principal camera, and multistereo frame is placed on its left
            % camera
            % - multicalib: 4x4xN where multicalib(:, :, ii) brings points
            % from camii to cam1
            n = length(stereos);
            obj.N = n;
            obj.Kls = zeros(3, 3, n);
            obj.Krs = zeros(3, 3, n);
            obj.Rs = zeros(3, 3, n);
            obj.ts = zeros(3, n);
            obj.stereosR = zeros(3, 3, n);
            obj.stereost = zeros(3, n);
            for ii = 1:n
                [Kl, Kr, Rstereo, tstereo] = obj.convertStereoParams(stereos{ii});
                obj.Kls(:, :, ii) = Kl;
                obj.Krs(:, :, ii) = Kr;
                obj.stereosR(:, :, ii) = Rstereo;
                obj.stereost(:, ii) = tstereo;
                obj.Rs(:, :, ii) = multicalib(1:3, 1:3, ii)';
                obj.ts(:, ii) = -multicalib(1:3, 1:3, ii)' * multicalib(1:3, 4, ii);
                obj.Pls(:, :, ii) = Kl * [obj.Rs(:, :, ii), obj.ts(:, ii)];
                obj.Prs(:, :, ii) = Kr * [Rstereo, tstereo] * [obj.Rs(:, :, ii), obj.ts(:, ii); zeros(1, 3), 1];
            end
        end
        
        function [Kl, Kr] = getKs(obj, cam_index)
            Kl = obj.Kls(:, :, cam_index);
            Kr = obj.Krs(:, :, cam_index);
        end
        
        function [R, t] = stereoLeft2Right(obj, cam_index)
            R = obj.stereosR(:, :, cam_index);
            t = obj.stereost(:, cam_index);
        end
        
        function [Rl, tl, Rr, tr] = Tprincipal2cam(obj, cam_index)
            % return transformation bringing points from principal camera
            % to the camera corresponding to cam_index
            Rl = obj.Rs(:, :, cam_index); 
            tl = obj.ts(:, cam_index);
            if nargout > 2
                Rr = obj.stereosR(:, :, cam_index) * Rl; 
                tr = obj.stereosR(:, :, cam_index) * tl + obj.stereost(:, cam_index);
            end
        end
        
        function [Rl, tl, Rr, tr] = Tcam2principal(obj, cam_index)
            % return transformation bringing points from the camera 
            % corresponding to cam_index to principal camera
            if nargout > 2
                [Rl_, tl_, Rr_, tr_] =  obj.Tprincipal2cam(cam_index);
                Rl = Rl_'; tl = - Rl_' * tl_;
                Rr = Rr_'; tr = - Rr_' * tr_;
            else
                [Rl_, tl_] =  obj.Tprincipal2cam(cam_index);
                Rl = Rl_'; tl = - Rl_' * tl_;
            end
        end
        
        function [Pl, Pr] = projMatrix(obj, cam_index)
            % return projection matrices in principal camera frame
            Pl = obj.Pls(:, :, cam_index);
            if nargout > 1
                Pr = obj.Prs(:, :, cam_index);
            end
        end
        
        function n = stereoPairsNumber(obj)
            n = obj.N;
        end
    end
    
    methods (Access = private)
        function [Kl, Kr, Rs, ts] = convertStereoParams(obj, s)
            Kl = s.CameraParameters1.IntrinsicMatrix';
            Kr = s.CameraParameters2.IntrinsicMatrix';
            Rs = s.RotationOfCamera2';
            ts = s.TranslationOfCamera2';
        end
    end
end

