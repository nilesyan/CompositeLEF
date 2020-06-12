function regenerate_event(LEFSystem, evheap, event_idx)

%     Regenerate an event in an event heap. If the event is currently impossible (e.g. a step
%     onto an occupied site) then the new event is not created, but the existing event is not
%     modified.
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


global ev_t14 ev_t13 ev_t6 ev_t7

% A step to the left or to the right.
if event_idx <= 14 * LEFSystem.N
    if event_idx <= 7 * LEFSystem.N
        leg_idx = event_idx;
        direction = -1;
    else
        leg_idx = event_idx - 7 * LEFSystem.N;
        direction = 1;
    end
    
    if LEFSystem.locs(leg_idx) > 0
        % Local velocity = velocity * permeability
        local_vel = (LEFSystem.perms(LEFSystem.locs(leg_idx) + (direction+1)/2)...
            * LEFSystem.vels(leg_idx + (direction*3.5+3.5) * LEFSystem.N));
        local_vel7 = (LEFSystem.perms(LEFSystem.locs(7) + (direction+1)/2)...
            * LEFSystem.vels(7));
        local_vel13 = (LEFSystem.perms(LEFSystem.locs(6) + (direction+1)/2)...
            * LEFSystem.vels(13));
        if local_vel > 0
            if (event_idx > 12*LEFSystem.N && event_idx < 13*LEFSystem.N+1 ...
                    && LEFSystem.lattice(LEFSystem.locs(leg_idx) + 1) > 0) ...
                    || (event_idx > 6*LEFSystem.N && event_idx < 7*LEFSystem.N+1 ...
                    && LEFSystem.lattice(LEFSystem.locs(leg_idx) - 1) == 6)
                evheap.remove_event(7);
                evheap.remove_event(13);
                ev_t13(1)=LEFSystem.time + exprnd(1/local_vel13/exp(LEFSystem.dG*(LEFSystem.f-1)));
                ev_t7(1)=LEFSystem.time + exprnd(1/local_vel7/exp(LEFSystem.dG*(LEFSystem.f-1)));
                evheap.add_event(17,ev_t7(1));
                evheap.add_event(18,ev_t13(1));
            elseif (event_idx > 4*LEFSystem.N && event_idx < 5*LEFSystem.N+1 ...
                    && LEFSystem.locs(leg_idx+1)==LEFSystem.locs(leg_idx)+2 && ...
                    LEFSystem.locs(leg_idx+2)==LEFSystem.locs(leg_idx)+2)
                
                local_vell = (LEFSystem.perms(LEFSystem.locs(6)) * LEFSystem.vels(6));
                ev_t6(1)=LEFSystem.time + exprnd(1/local_vell/exp(LEFSystem.dG*(LEFSystem.f)));
                evheap.add_event(21,ev_t6(1));
                if LEFSystem.lattice(LEFSystem.locs(leg_idx) + direction) <= 0
                    ev_t=LEFSystem.time + exprnd(1/local_vel);
                    evheap.add_event(event_idx,ev_t);
                end
            else
                if LEFSystem.lattice(LEFSystem.locs(leg_idx) + direction) <= 0
                    ev_t=LEFSystem.time + exprnd(1/local_vel);
                    if event_idx > 13*LEFSystem.N && event_idx < 14*LEFSystem.N+1
                        ev_t14(event_idx-13*LEFSystem.N)=ev_t;
                    end
                    if event_idx > 5*LEFSystem.N && event_idx < 6*LEFSystem.N+1
                        ev_t6(event_idx-5*LEFSystem.N)=ev_t;
                    end
                    evheap.add_event(event_idx,ev_t);
                end
            end
        end
    end
    
    % Passive unbinding.
elseif event_idx > 14 * LEFSystem.N && event_idx <= 15 * LEFSystem.N
    loop_idx = event_idx - 14 * LEFSystem.N;
    if LEFSystem.locs(loop_idx) > 0 && LEFSystem.locs(loop_idx+LEFSystem.N) > 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.lifespans(loop_idx)))
    end
    
    % Rebinding from the solution to a random site.
elseif event_idx > 15 * LEFSystem.N && event_idx <= 16 * LEFSystem.N
    loop_idx = event_idx - 15 * LEFSystem.N;
    if LEFSystem.locs(loop_idx) <= 0 && LEFSystem.locs(loop_idx+LEFSystem.N) <= 0
        evheap.add_event(event_idx,LEFSystem.time + exprnd(LEFSystem.rebinding_times(loop_idx)))
    end
    
elseif event_idx > 18 * LEFSystem.N && event_idx <= 19 * LEFSystem.N
    if LEFSystem.locs(6)==LEFSystem.locs(7)
        evheap.remove_event(13);
        evheap.remove_event(7);
        evheap.remove_event(14);
        evheap.remove_event(6);
        local_velr = (LEFSystem.perms(LEFSystem.locs(7) + 1)...
            * LEFSystem.vels(14));
        local_vell = (LEFSystem.perms(LEFSystem.locs(6) )...
            * LEFSystem.vels(6));
        ev_t14(1)=LEFSystem.time + exprnd(1/local_velr/exp(LEFSystem.dG*(LEFSystem.f)));
        ev_t6(1)=LEFSystem.time + exprnd(1/local_vell/exp(LEFSystem.dG*(LEFSystem.f)));
        evheap.add_event(20,ev_t14(1));
        evheap.add_event(21,ev_t6(1));
    end
end

end

