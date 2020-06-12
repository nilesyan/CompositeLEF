classdef Event_t < handle
    properties
        time
        event_idx
    end
    methods
        function obj = Event_t(time, event_idx)
            obj.time = time;
            obj.event_idx = event_idx;
        end    
    end
end