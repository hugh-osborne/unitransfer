timestep = 0.1;

g_nap = 0.25; %mS - %250.0; %nS
theta_m = -47.1; %mV
sig_m = -3.1; %mV
theta_h = -59; %mV
sig_h = 8; %mV
tau_h = 1200; %ms 
E_na = 55; %mV
C = 1; %uF - %1000; %pF
% Depending on how much noise we expect, E_l should be reduced to the point
% where there is still no activity - that is, the stationary point comes before the v
% nullcline peak. (higher I => lower v nullcline peak - stationary point moves to the right, 
% lower E_l => higher nullcline peak - stationary point moves to the left)
% ***Note for each -1 in E_l, I must increase by g_l to keep the same stationary point.*** 
g_l = 0.1; %mS
E_l = -64.0; %mV
I = 0.0; %
I_h = 0; %

start_point_x = -59.5;

start_point_y = ((I - g_l.*(start_point_x-E_l)).*(1 + exp((start_point_x-theta_m)/sig_m)))./(g_nap.*(start_point_x-E_na));

tspan = 0:timestep:3;
[t1,s1] = ode23s(@izshikevich, tspan, [start_point_x start_point_y]);

x1 = s1(:,1);
y1 = s1(:,2);

centre_point_x = -52.85;
time_end = min(find(x1>centre_point_x))*timestep;
time_55 = min(find(y1>0.55))*timestep;

plot(x1,y1,'k');
hold on;

nvs = -80:0.1:0;
    
v_nulls = ((I - g_l.*(nvs-E_l)).*(1 + exp((nvs-theta_m)/sig_m)))./(g_nap.*(nvs-E_na));
h_nulls = 1 ./ (1+exp((nvs-theta_h)./sig_h));

plot(nvs,v_nulls,'r');
hold on;

plot(nvs,h_nulls,'b');
hold on;

I = 0.25; %

v_nulls = ((I - g_l.*(nvs-E_l)).*(1 + exp((nvs-theta_m)/sig_m)))./(g_nap.*(nvs-E_na));
h_nulls = 1 ./ (1+exp((nvs-theta_h)./sig_h));

plot(nvs,v_nulls,'g');
hold on;

plot(nvs,h_nulls,'c');
hold on;

stat_min_v = -59.05;
stat_max_v = -58.95;
stat_min_w = 0.56;
stat_max_w = 0.54;

axis fill
xlim([-70 -40]);
ylim([0 1.2]);
    
title('Rybak P-Na Neuron');
ylabel('h');
xlabel('V');