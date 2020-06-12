function sta = do_event(LEFSystem, evheap, event_idx)

%     Apply an event from a heap on the system and then regenerate it.
%     If the event is currently impossible (e.g. a step onto an occupied site),
%     it is not applied, however, no warning is raised.
%
%     Also, partially checks the system for consistency. Returns 0 if the system
%     is not consistent (a very bad sign), otherwise returns 1 if the event was a step
%     and 2 if the event was rebinding.
%
%     Possible events:
%     1 to 7N : a step to the left
%     7N+1 to 14N : a step to the right
%     14N+1 to 15N : passive unbinding
%     15N+1 to 16N : rebinding to a randomly chosen site
%     16N+1 to 17N : #7 step left onto #6
%     17N+1 to 18N : #6 step right onto #7
%     18N+1 to 19N : when #6 & #7 are together
%     19N+1 to 20N : #7 leaves and step right
%     20N+1 to 21N : #6 leaves and step left


if event_idx <= 14 * LEFSystem.N
    % Take a step
    if event_idx <= 7 * LEFSystem.N
        leg_idx = event_idx;
        direction = -1;
    else
        leg_idx = event_idx - 7 * LEFSystem.N;
        direction = 1;
    end
    
    prev_pos = LEFSystem.locs(leg_idx);
    % check if the loop was attached to the chromatin
    sta = 1;
    if prev_pos > 0
        % make a step only if there is no boundary and the new position is unoccupied
        if LEFSystem.perms(prev_pos + (direction + 1) / 2) > 0
            if LEFSystem.lattice(prev_pos+direction) < 0
                sta = sta * LEFSystem.make_step(leg_idx, direction);
                % regenerate events for the previous and the new neighbors
                if (event_idx>12*LEFSystem.N && event_idx < 13*LEFSystem.N+1 && LEFSystem.lattice(prev_pos+1)>0)...
                        || (event_idx>6*LEFSystem.N && event_idx < 7*LEFSystem.N+1 && LEFSystem.lattice(prev_pos-1)>0)
                    regenerate_neighbours(LEFSystem, evheap, prev_pos);
                elseif event_idx>4*LEFSystem.N && event_idx < 5*LEFSystem.N+1 &&...
                        LEFSystem.locs(leg_idx+1)==LEFSystem.locs(leg_idx)+2 && ...
                        LEFSystem.locs(leg_idx+2)==LEFSystem.locs(leg_idx)+2
                    regenerate_neighbours(LEFSystem, evheap, prev_pos+direction);
                    regenerate_event(LEFSystem, evheap, 12)
                    
                elseif (event_idx > 11*LEFSystem.N && event_idx < 12*LEFSystem.N+1 ...
                        && LEFSystem.locs(leg_idx+1)==LEFSystem.locs(leg_idx)+1 && ...
                        LEFSystem.locs(leg_idx+2)==LEFSystem.locs(leg_idx)+1)
                    regenerate_neighbours(LEFSystem, evheap, prev_pos);
                    regenerate_event(LEFSystem, evheap, 5)
                    
                else
                    regenerate_neighbours(LEFSystem, evheap, prev_pos);
                    regenerate_neighbours(LEFSystem, evheap, prev_pos+direction);
                    
                end
                
            end
        end
        
        % regenerate the performed event
        regenerate_event(LEFSystem, evheap, event_idx);
    end
    
elseif event_idx > 14 * LEFSystem.N && event_idx <= 15 * LEFSystem.N
    % unbinding
    loop_idx = event_idx - 14 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was attached to the chromatin
    if LEFSystem.locs(loop_idx) <= 0 || LEFSystem.locs(loop_idx+LEFSystem.N) <= 0
        sta = 0;
    end
    
    % save previous positions, but don't update neighbours until the loop
    % has moved
    prev_pos1 = LEFSystem.locs(loop_idx);
    prev_pos2 = LEFSystem.locs(loop_idx + LEFSystem.N);
    
    sta = sta * LEFSystem.move_leg(loop_idx, -1);
    sta = sta * LEFSystem.move_leg(loop_idx+LEFSystem.N, -1);
    
    % regenerate events for the loop itself and for its previous neighbours
    regenerate_all_loop_events(LEFSystem, evheap, loop_idx);
    
    % update the neighbours after the loop has moved
    regenerate_neighbours(LEFSystem, evheap, prev_pos1);
    regenerate_neighbours(LEFSystem, evheap, prev_pos2);
    
elseif event_idx > 15 * LEFSystem.N && event_idx <= 16 * LEFSystem.N
    loop_idx = event_idx - 15 * LEFSystem.N;
    
    sta = 2;
    % check if the loop was not attached to the chromatin
    if LEFSystem.locs(loop_idx) > 0 || LEFSystem.locs(loop_idx+LEFSystem.N) > 0
        sta = 0;
    end
    
    % find a new position for the LEM (a brute force method, can be
    % improved)
    while 1
        new_pos = randi([1 LEFSystem.L-1],1,1);
        if LEFSystem.lattice(new_pos) <= 0 && LEFSystem.lattice(new_pos+1) <= 0 && LEFSystem.perms(new_pos+1) > 0
            break
        end
    end
    % rebind the loop
    sta = sta * LEFSystem.move_leg(loop_idx, new_pos);
    sta = sta * LEFSystem.move_leg(loop_idx+LEFSystem.N, new_pos+1);
    
    % regenerate events for the loop itself and for its new neighbours
    regenerate_all_loop_events(LEFSystem, evheap, loop_idx);
    regenerate_neighbours(LEFSystem, evheap, new_pos);
    regenerate_neighbours(LEFSystem, evheap, new_pos + 1);
    
elseif event_idx > 16 * LEFSystem.N && event_idx <= 17 * LEFSystem.N
    %loop_idx = event_idx - 16 * LEFSystem.N;
    leg_idx=7;
    prev_pos = LEFSystem.locs(leg_idx);
    if  LEFSystem.locs(leg_idx-1)==prev_pos-1
        LEFSystem.lattice(prev_pos)=-1;
        LEFSystem.locs(leg_idx)=LEFSystem.locs(leg_idx)-1;
        evheap.remove_event(18);
        regenerate_event(LEFSystem, evheap,19);
    end
    %regenerate_event(LEFSystem, evheap,20);
    sta=1;
    
elseif event_idx > 17 * LEFSystem.N && event_idx <= 18 * LEFSystem.N
    leg_idx=6;
    prev_pos = LEFSystem.locs(leg_idx);
    if LEFSystem.locs(leg_idx+1)==prev_pos+1
        LEFSystem.locs(leg_idx)=LEFSystem.locs(leg_idx)+1;
        LEFSystem.lattice(prev_pos)=-1;
        evheap.remove_event(17);
        regenerate_event(LEFSystem, evheap,19);
        if prev_pos > 1 && LEFSystem.lattice(prev_pos-1) > 0
            regenerate_event(LEFSystem, evheap, LEFSystem.lattice(prev_pos-1) + 7 * LEFSystem.N)
        end
    end
    sta=1;
    
elseif event_idx > 19 * LEFSystem.N && event_idx <= 20 * LEFSystem.N
    leg_idx=7;
    prev_pos = LEFSystem.locs(leg_idx);
    if LEFSystem.locs(leg_idx)==LEFSystem.locs(leg_idx-1)
        LEFSystem.make_step(leg_idx, 1);
        LEFSystem.lattice(prev_pos)=leg_idx-1;
        regenerate_event(LEFSystem, evheap,6);
        regenerate_event(LEFSystem, evheap,14);
        regenerate_event(LEFSystem, evheap,7);
        regenerate_event(LEFSystem, evheap,13);
        evheap.remove_event(21);
    end
    sta=1;
    
elseif event_idx > 20 * LEFSystem.N && event_idx <= 21 * LEFSystem.N
    leg_idx=6;
    prev_pos = LEFSystem.locs(leg_idx);
    if LEFSystem.locs(leg_idx)==LEFSystem.locs(leg_idx+1)
        if LEFSystem.lattice(prev_pos-1)<0
            LEFSystem.make_step(leg_idx, -1);
            LEFSystem.lattice(prev_pos)=leg_idx+1;
            regenerate_event(LEFSystem, evheap,6);
            regenerate_event(LEFSystem, evheap,14);
            regenerate_event(LEFSystem, evheap,7);
            regenerate_event(LEFSystem, evheap,13);
            evheap.remove_event(20);
            if LEFSystem.locs(leg_idx-1)==prev_pos-2
                regenerate_event(LEFSystem, evheap,12);
            end
        end
    end
    sta=1;
    
else
    disp(strcat('event_idx assumed a forbidden value :', num2str(event_idx)))
    sta = 0;
end
end