function regenerate_neighbours(LEFSystem, evheap, pos)

%     Regenerate the motion events for the adjacent loop legs.
%     Use to unblock the previous neighbors and block the new ones.

% regenerate left step for the neighbour on the right
if pos < LEFSystem.L && LEFSystem.lattice(pos+1) > 0
    regenerate_event(LEFSystem, evheap, LEFSystem.lattice(pos+1))
end

% regenerate right step for the neighbour on the left
if pos > 1 && LEFSystem.lattice(pos-1) > 0
    regenerate_event(LEFSystem, evheap, LEFSystem.lattice(pos-1) + 7 * LEFSystem.N)
end

end