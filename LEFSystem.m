classdef LEFSystem < handle
    properties
        L
        N
        time
        vels
        lifespans
        rebinding_times
        perms
        lattice
        locs
        dG
        f
    end
    methods
        function obj = LEFSystem(L,N,vels,lifespans,rebinding_times,init_locs,perms,BE,BE_perms,dG,f)
            
            % the occupancy of each site of the system. If -1, the site is unoccupied.
            obj.L = L;
            obj.N = N;
            obj.lattice = -1 * ones(L,1);
            obj.locs = -1 * ones(7*N,1);
            obj.vels = vels;
            obj.lifespans = lifespans;
            obj.rebinding_times = rebinding_times;
            obj.dG = dG;
            obj.f = f;
            
            if isempty(BE)
                if isempty(perms)
                    obj.perms = ones(L+1,1);
                else
                    obj.perms=perms;
                end
            else
                if isempty(perms)
                    obj.perms = ones(L+1,1);
                else
                    obj.perms=perms;
                end
                for i = 1:length(BE)
                    obj.perms(BE(i)+1) = BE_perms;
                end
            end
            
            obj.perms(1)=0;
            obj.perms(end)=0;
            
            % Initialize non-random loops
            for i = 1:7*N
                % Check if the loop is preinitialized.
                if init_locs(i) < 0
                    continue
                end
                % Populate a site.
                obj.locs(i) = init_locs(i);
                obj.lattice(obj.locs(i)) = i;
            end
        end
        
        function r = make_step(obj,leg_idx,direction)
            % The variable `direction` can only take values +1 or -1.
            new_pos = obj.locs(leg_idx) + direction;
            r = obj.move_leg(leg_idx, new_pos);
        end
        
        function r = move_leg(obj,leg_idx,new_pos)
            if new_pos > 0 && obj.lattice(new_pos) > 0
                r = 0;
                return
            end
            prev_pos = obj.locs(leg_idx);
            obj.locs(leg_idx) = new_pos;
            
            if prev_pos > 0
                if obj.lattice(prev_pos) <= 0
                    r = 0;
                    return
                end
                obj.lattice(prev_pos) = -1;
            end
            if new_pos > 0
                obj.lattice(new_pos) = leg_idx;
            end
            r = 1;
        end
        
        function okay = check_system(obj)
            okay = 1;
            for i = 1:obj.N
                if obj.locs(i) == obj.locs(i+obj.N)
                    disp(strcat('loop ' , num2str(i), 'has both legs at ', num2str(obj.locs(i))))
                    okay = 0;
                end
                if obj.locs(i) > obj.L
                    disp(strcat('leg ', num2str(i), 'is located outside of the system: ', num2str(obj.locs(i))))
                    okay = 0;
                end
                if obj.locs(i+obj.N) > obj.L
                    disp(strcat('leg ', num2str(i+obj.N), 'is located outside of the system: ', num2str(obj.locs(i+obj.N))))
                    okay = 0;
                end
                if (obj.locs(i) <= 0 && obj.locs(i+obj.N) > 0 ) || (obj.locs(i) > 0 && obj.locs(i+obj.N) <= 0 )
                    disp(strcat('the legs of the loop', num2str(i), 'are inconsistent: ', num2str((obj.locs(i))), num2str(obj.locs(i+obj.N))))
                    okay = 0;
                end
            end
        end
    end
end