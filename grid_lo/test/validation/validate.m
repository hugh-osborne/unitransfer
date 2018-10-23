axis fill
xlim([-90 -40]);
ylim([-200 200]);

title('Grid');
ylabel('v');
xlabel('t');

timestep = 1;
tspan = 0:timestep:1000;
last_time = 0;
last_h = 0.0;
global spike_count;
spike_count = 0;
reset = odeset('Events', @freset);

[t,s,te,ye,ie] = ode23s(@adex, tspan, [-72 0.0], reset);
vs = s(:,1);
last_h = s(end,2);
plot(t,vs,'r');
last_time = te;
hold on;
tspan = te:timestep:1000.0;

while (last_time < 1000)
    [t,s,te,ye,ie] = ode23s(@adex, tspan, [-58 last_h], reset);
    vs = s(:,1);
    last_h = s(end,2);
    last_time = te;
    spike_count = spike_count + 1;
    plot(t,vs,'r');
    hold on;
    tspan = te:timestep:1000.0;
end

spike_count


function [value, isterminal, direction] = freset(t, vs)
    value      = (vs(1) > -10)-0.5;
    isterminal = 1;
    direction  = 0;
end

function r = adex(t,y)
    C = 200;
    g_l = 10;
    E_l = -70.0;
    v_t = -50.0;
    tau = 2.0;
    alpha = 2.0;
    tau_w = 30.0;

    v = y(1);
    w = y(2);
    
    input = poissrnd(5000,1,1);
    add = input(1) * 0.1
    
    v_prime = (-g_l*(v - E_l) + g_l*tau*exp((v - v_t)/tau) - w + add) / C;
    w_prime = (alpha*(v-E_l) - w) / tau_w;

    r = [v_prime; w_prime];
end