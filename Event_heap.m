classdef Event_heap < handle
    properties
        heap
        entry_finder
    end
    methods
        function obj = Event_heap()
            obj.heap = containers.Map('KeyType','double','ValueType','any');
            obj.entry_finder = containers.Map('KeyType','int64','ValueType','any');
        end
        
        function add_event(obj,event_idx,time)
            % Add a new event or update the time of an existing event.
            if isKey(obj.entry_finder,event_idx)
                obj.remove_event(event_idx)
            end
            entry = Event_t(time, event_idx);
            obj.entry_finder(event_idx) = entry;
            obj.heap(entry.event_idx) = entry;
        end
        
        function remove_event(obj,event_idx)
            % Mark an existing event as REMOVED.
            if isKey(obj.entry_finder,event_idx)
                entry = obj.entry_finder(event_idx);
                remove(obj.entry_finder,event_idx);
                remove(obj.heap,event_idx);
                %entry.event_idx = -1;
            end
        end
        
        function outevent = pop_event(obj)
            % Remove and return the closest event.
            while size(obj.heap,1)
                values = obj.heap.values;
                for i = 1:size(obj.heap,1)
                    sorttime(i)=values{i}.time;
                    sortevent(i)=values{i}.event_idx;
                    if sortevent(i)<0
                        disp(009)
                    end
                end
                [~,sortid] = min(sorttime);
                entry = obj.heap(sortevent(sortid));
                remove(obj.heap,sortevent(sortid));
                if entry.event_idx ~= -1
                    remove(obj.entry_finder,entry.event_idx);
                    outevent = entry;
                    return
                end
            end
            outevent = Event_t(0, 0);
        end
        
    end
end