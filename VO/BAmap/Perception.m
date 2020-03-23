classdef Perception
    %PERCEPTION An image perception of a 3D point associated with timestamp and camera index
    
    properties
        P2D     % image point
        CamIdx  % index of the perceiving camera
        Ts      % timestamp at which the point is perceived
    end
    
    methods
        function obj = Perception(p2D, ts, camIdx)
            %PERCEPTION Construct an instance of this class
            obj.P2D = p2D;
            obj.Ts = ts;
            obj.CamIdx = camIdx;
        end
    end
end

