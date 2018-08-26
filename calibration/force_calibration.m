
%% from Chaigne and Askenfelt II
SR = 44.1e3;
tr = [1.5e-3, 2.5e-3]; % milliseconds, time range
fr = [10 40];    % force range, Newtons.
vs = [1 127];    % range of velocities

%% Creating force profile, the softer the force, the longer the time
v = 1;
newtons = (fr(2)-fr(1))*(v-1)/(vs(2)-vs(1)) + fr(1)
dur = floor((((tr(2)-tr(1))*(newtons - fr(1))/(fr(2) - fr(1)))+tr(1))*SR);                 % duration in samples
tn = 0:1:(dur);                 % time vector in samples

rc = (1-cos(2*pi*(tn/dur)))*.5; %% raised cosine

plot(rc)

%%% scaling one range to another
%% f(x) = (b-a)(x-min)/(max-min) + a
