function regenerate_all_loop_events(LEFSystem, evheap, loop_idx)

%     Regenerate all possible events for a loop. Includes the four possible motions
%     and passive unbinding.

    regenerate_event(LEFSystem, evheap, loop_idx)
    regenerate_event(LEFSystem, evheap, loop_idx + LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 2 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 3 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 4 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 5 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 6 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 7 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 8 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 9 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 10 * LEFSystem.N)  
    regenerate_event(LEFSystem, evheap, loop_idx + 11 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 12 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 13 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 14 * LEFSystem.N)
    regenerate_event(LEFSystem, evheap, loop_idx + 15 * LEFSystem.N)    
end