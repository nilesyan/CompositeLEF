clear
clc
%% Parameters
rng(1)
L = 500; %lattice length
N = 1; %number of LEF pairs
km = 1e-2; %remodeler back-stepping rate
kp = 0.2; %remodeler forward-stepping rate
beta = 10; %nucleosome unbinding rate
alpha = 133; %nucleosome binding rate
m = 1; %SMC stepping rate (here assumes symmetric rate)
f = 1;
deltag = log(alpha/beta);
dG = 8;
R_EXTEND = [beta km alpha beta m km alpha];
R_SHRINK = [alpha kp beta alpha m kp beta];
%R_EXTEND & R_SHRINK: The last 4 components form a composite LEF foot. The
%first 3 form an incomplete LEF foot and don't interfere with the last 4.
%(fur further study purpose).
%[junction remodeler junction junction SMC remodeler junction]
R_OFF = 1e-200; %fall-off rate. An extremely small rate is used to avoid fall off.
REBINDING_TIME = 0;
INIT_SITES = [1]; %Initial positions are assigned below.
ACTIVATION_TIMES = zeros(N,1);
T_MAX = 5e0;
N_SNAPSHOTS = 1e5;
BE = [];
BE_perms = 0;
PERMS = [];
verbose = 1;

%% Initiation
global ev_t14 ev_t13 ev_t6 ev_t7
ev_t14= -1*ones(N);
ev_t13= -1*ones(N);
ev_t6= -1*ones(N);
ev_t7= -1*ones(N);

VELS = zeros(14*N,1);
LIFESPANS = zeros(N,1);
REBINDING_TIMES = zeros(N,1);

for i = 1:N
    if length(R_EXTEND) > 1
        VELS(i) = R_EXTEND(1);
        VELS(i+N) = R_EXTEND(2);
        VELS(i+2*N) = R_EXTEND(3);
        VELS(i+3*N) = R_EXTEND(4);
        VELS(i+4*N) = R_EXTEND(5);
        VELS(i+5*N) = R_EXTEND(6);
        VELS(i+6*N) = R_EXTEND(7);
    else
        VELS(i) = R_EXTEND;
        VELS(i+3*N) = R_EXTEND;
    end
    
    if length(R_SHRINK) > 1
        VELS(i+7*N) = R_SHRINK(1);
        VELS(i+8*N) = R_SHRINK(2);
        VELS(i+9*N) = R_SHRINK(3);
        VELS(i+10*N) = R_SHRINK(4);
        VELS(i+11*N) = R_SHRINK(5);
        VELS(i+12*N) = R_SHRINK(6);
        VELS(i+13*N) = R_SHRINK(7);
    else
        VELS(i+N) = R_SHRINK;
        VELS(i+2*N) = R_SHRINK;
    end
    
    if length(R_OFF) > 1
        LIFESPANS(i) = 1/R_OFF(i);
    else
        LIFESPANS(i) = 1/R_OFF;
    end
    
    if length(REBINDING_TIME) > 1
        REBINDING_TIMES(i) = REBINDING_TIME(i);
    else
        REBINDING_TIMES(i) = REBINDING_TIME;
    end
end

INIT_LOCS = -1 * ones(7*N,1);

for i = 1:N
    if INIT_SITES(i)>0
        INIT_LOCS(i) = INIT_SITES(i);
        INIT_LOCS(i+N) = INIT_SITES(i)+1;
        INIT_LOCS(i+2*N) = INIT_SITES(i)+2;
        % The first 3 form an incomplete LEF foot that starts at 1.
        INIT_LOCS(i+3*N) = INIT_SITES(i)+203;
        INIT_LOCS(i+4*N) = INIT_SITES(i)+204;
        INIT_LOCS(i+5*N) = INIT_SITES(i)+205;
        INIT_LOCS(i+6*N) = INIT_SITES(i)+206;
        % The last 4 form a complete LEF foot that starts at a far distance
        % from the other foot. So the two feet will not interfere.
    end
end

if ~isempty(PERMS) && length(PERMS) ~= L+1
    msgID = 'PERMS:Length';
    msg = 'The length of the provided array of permeabilities should be L+1.';
    Exception = MException(msgID,msg);
    throw(Exception)
end

LEFSYSTEM = LEFSystem(L, N, VELS, LIFESPANS, REBINDING_TIMES, INIT_LOCS,...
    PERMS, BE, BE_perms, dG, f);
LEFSYSTEM.time=0;

sites_traj1 = zeros(N_SNAPSHOTS, N);
sites_traj2 = zeros(N_SNAPSHOTS, N);
sites_traj3 = zeros(N_SNAPSHOTS, N);
sites_traj4 = zeros(N_SNAPSHOTS, N);
sites_traj5 = zeros(N_SNAPSHOTS, N);
sites_traj6 = zeros(N_SNAPSHOTS, N);
sites_traj7 = zeros(N_SNAPSHOTS, N);
ts_traj = zeros(N_SNAPSHOTS, 1);

prev_snapshot_t = 0;
snapshot_idx = 1;

evheap = Event_heap();

for i = 1:LEFSYSTEM.N
    if INIT_LOCS(i) == -1 && INIT_LOCS(i) == -1
        evheap.add_event(i + 5 * LEFSYSTEM.N, ACTIVATION_TIMES(i));
    else
        regenerate_all_loop_events(LEFSYSTEM, evheap, i)
        regenerate_neighbours(LEFSYSTEM, evheap, INIT_LOCS(i))
        regenerate_neighbours(LEFSYSTEM, evheap, INIT_LOCS(i+LEFSYSTEM.N))
    end
end

%% Iteration
while snapshot_idx <= N_SNAPSHOTS
    LEFEvent = evheap.pop_event();
    LEFSYSTEM.time = LEFEvent.time;
    event_idx = LEFEvent.event_idx;
    LEFStatus = do_event(LEFSYSTEM, evheap, event_idx);
    
    if LEFStatus == 0
        disp('an assertion failed somewhere')
        return
    end
    
    if LEFSYSTEM.time > prev_snapshot_t + T_MAX / N_SNAPSHOTS
        prev_snapshot_t = LEFSYSTEM.time;
        sites_traj1(snapshot_idx,1:N) = LEFSYSTEM.locs(1:N);
        sites_traj2(snapshot_idx,1:N) = LEFSYSTEM.locs(N+1:2*N);
        sites_traj3(snapshot_idx,1:N) = LEFSYSTEM.locs(2*N+1:3*N);
        sites_traj4(snapshot_idx,1:N) = LEFSYSTEM.locs(3*N+1:4*N);
        sites_traj5(snapshot_idx,1:N) = LEFSYSTEM.locs(4*N+1:5*N);
        sites_traj6(snapshot_idx,1:N) = LEFSYSTEM.locs(5*N+1:6*N);
        sites_traj7(snapshot_idx,1:N) = LEFSYSTEM.locs(6*N+1:end);
        ts_traj(snapshot_idx) = LEFSYSTEM.time;
        snapshot_idx = snapshot_idx + 1;
                if verbose && mod(snapshot_idx,1e4) == 0
                    disp([snapshot_idx/N_SNAPSHOTS])
                end
    end
end
v_rem= (sites_traj6(end)-sites_traj6(1))/(ts_traj(end)-ts_traj(1)); %remodeler velocity
v_SMC = (sites_traj5(end)-sites_traj5(1))/(ts_traj(end)-ts_traj(1)); %SMC velocity
D_SMC = var(sites_traj5)/2/ts_traj(end); %remodeler diffusivity
D_rem = var(sites_traj6)/2/ts_traj(end); %SMC diffusivity

%% Plot of the complete LEF foot
ist = 200;
figure;
plot(ts_traj,sites_traj5-ist,'r','linewidth',2);
hold on;
plot(ts_traj,sites_traj6-ist,'b','linewidth',2);
plot(ts_traj,sites_traj7-ist,'color',[0.8 0.8 0.8]);
plot(ts_traj,sites_traj4-ist,'color',[0.8 0.8 0.8]);
xlabel('Time')
ylabel('Position')
legend('SMC Protein','Polymerase','Nucleosome Junctions','Location','northwest')
set(gca,'FontSize', 18)

%% Plot of the incomplete LEF foot
ist = INIT_SITES;
figure;plot(ts_traj,sites_traj2-ist,'r','linewidth',2);
hold on;
plot(ts_traj,sites_traj1-ist,'color',[0.8 0.8 0.8]);
plot(ts_traj,sites_traj3-ist,'color',[0.8 0.8 0.8]);
xlabel('Time')
ylabel('Position')
legend('SMC Protein','Nucleosome Juntion','Location','northwest')
set(gca,'FontSize', 18)

%% Histogram of separation
Sep = sites_traj6-sites_traj5; %separation between remodeler and SMC
figure;
hist(Sep)