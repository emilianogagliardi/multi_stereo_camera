classdef Landmark
    % LADMARK a 3D point with a list of 2D perceptions
    
    properties
        P3D
        Perceptions
    end
    
    methods
        function obj = Landmark(p3D, p2D, ts, camIdx)
            %LADMARK Construct an instance of this class
            perception = Perception(p2D, ts, camIdx);
            obj.P3D = p3D;
            obj.Perceptions = [perception];
        end
        
        function obj = AddPerception(obj, p2D, ts, camIdx)
            perception = Perception(p2D, ts, camIdx);
            obj.Perceptions = [obj.Perceptions, perception];
        end
        
        function [percs, obj] = GetPerceptions(obj, from, to)
            percs_times = vertcat(obj.Perceptions.Ts);
            idx = percs_times >= from & percs_times <= to;
            percs = obj.Perceptions(idx);
        end
    end
end

