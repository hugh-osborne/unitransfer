global C;
global g_na;
global v_na;
global g_k;
global v_k;
global g_l;
global v_l;
global I;
    
C = 1;
g_na = 120;
v_na = 115;
g_k = 36;
v_k = -12;
g_l = 0.3;
v_l = -54;

v_0 = -64;
h_0 = h_inf(-64);
m_0 = m_inf(-64);
n_0 = n_inf(-64);

% Input current
I = 0;

% v_rest is the resting potential of this neuron
% User must find the resting potential :

global approx_a;
global approx_b;
global alpha;
global v_rest;

v_rest = -64; %-45

% alpha is the angle of the slope of the curve (n_inf, h_inf) at v_rest
% User must find the angle:

alpha = -1.38; %-0.785

approx_a = -tan(alpha);
approx_b = (approx_a * n_inf(v_rest)) + h_inf(v_rest);

% vs = -100:0.01:10;
% hs = h_inf(vs);
% ns = n_inf(vs);
% plot(ns,hs,'k');
% hold on;
% 
% zs = 0:0.01:1;
% hs = h_inf(v_rest) + zs*sin(alpha); 
% ns = n_inf(v_rest) + zs*cos(alpha); 
% plot(ns,hs,'r');


global z_0;
z_0 = (h_0 - h_inf(v_rest)) / sin(alpha);

%plot taus
% vspan = -100:0.1:100;
% ms = tau_m(vspan);
% hs = tau_h(vspan);
% ns = tau_n(vspan);
% plot(vspan, ms, 'r');
% hold on;
% plot(vspan, hs, 'g');
% plot(vspan, ns, 'b');

%plot infs
% vspan = -100:0.1:100;
% ms = m_inf(vspan);
% hs = h_inf(vspan);
% ns = n_inf(vspan);
% plot(vspan, ms, 'r');
% hold on;
% plot(vspan, hs, 'g');
% plot(vspan, ns, 'b');

% plot the full HH model
tspan = 0:0.01:50.0;
[t,s] = ode23s(@HHv_full, tspan, [v_0 m_0 h_0 n_0]);
vs = s(:,1);
plot(t,vs,'k');
hold on;
% % 3D using m_inf
% [t,s] = ode23s(@HHv_m_inf, tspan, [v_0 h_0 n_0]);
% vs = s(:,1);
% plot(t,vs,'b');
% % 3D using n h approximation
% [t,s] = ode23s(@HHv_nh, tspan, [v_0 m_0 h_0]);
% vs = s(:,1);
% plot(t,vs,'r');
% 2D using nh approximation and m_inf
% [t,s] = ode23s(@HHv_2D_nh_m_inf, tspan, [v_0 h_0]);
% vs = s(:,1);
% plot(t,vs,'b');
% 1D using nh approximation, m_inf anf h_inf
% tspan = 0:0.01:50.0;
% [t,s] = ode23s(@HHv_1D_nh_m_inf, tspan, [v_0]);
% vs = s(:,1);
% plot(t,vs,'g');

%plot 2D approximation cycle
% tspan = 0:0.01:50.0;
% [t,s] = ode23s(@HHv_2D_nh_m_inf, tspan, [v_0 h_0]);
% vs = s(:,1);
% hs = s(:,2);
% plot(vs,hs,'b');
% hold on;

%plot 1D approximation cycle
% tspan = 0:0.01:50.0;
% [t,s] = ode23s(@HHv_1D_nh_m_inf, tspan, [v_0]);
% vs = s(:,1);
% plot(vs,hs,'r');

%plot m
% tspan = 0:0.01:50.0;
% [t,s] = ode23s(@HHv_full, tspan, [v_0 m_0 h_0 n_0]);
% ms = s(:,2);
% plot(t,ms,'k');
% hold on;
% [t,s] = ode23s(@HHv_m_inf, tspan, [v_0 m_0 h_0 n_0]);
% m_infs = s(:,2);
% plot(t,m_infs,'r');


%plot h
% tspan = 0:0.01:50.0;
% [t,s] = ode23s(@HHv_full, tspan, [v_0 m_0 h_0 n_0]);
% hs = s(:,3);
% plot(t,hs,'k');
% hold on;
% [t,s] = ode23s(@HHv_h_inf, tspan, [v_0 m_0 h_0 n_0]);
% hs = s(:,3);
% plot(t,hs,'r');

%plot n and h
% use this to work out the linear approximation to h and n
% tspan = 0:0.01:100.0;
% [t,s] = ode23s(@HHv_full, tspan, [v_0 m_0 h_0 n_0]);
% hs = s(:,3);
% ns = s(:,4);
% plot(ns(2000:end),hs(2000:end),'b');
% hold on;

% with default values
% h/n step is (0.491,0.283) to (0.73,0.057)
% slope is (0.057 - 0.283) / (0.73 - 0.491) = -0.945
% h intercept is 0.283 + (0.491*0.945) = 0.747
% -> h = 0.747 - 0.945n
% -> n = 0.79 - 1.06h

%confirm the approximation
% [t,s] = ode23s(@HHv_nh, tspan, [v_0 m_0 h_0]);
% hs = s(:,3);
% ns = 0.79 - (1.06*hs);
% plot(ns(2000:end),hs(2000:end),'r');

% plot the 1D model
% tspan = 0:0.01:100.0;
% for i = 0:5
%     reset = odeset('Events', @HHv_reset);
%     [t,s,te,ye,ie] = ode23s(@HHv_1D, tspan, [v_0], reset);
%     vs = s(:,1);
%     plot(t,vs,'k');
%     hold on;
%     v_0 = -80;
%     h_0 = h_inf(-80);
%     m_0 = m_inf(-80);
%     n_0 = n_inf(-80);
%     tspan = te+5:0.01:100.0;
% end

% Generate 1D mesh

% global timestep;
% global revId;
% global outId;
% global strip;
% global formatSpec;
% 
% timestep = 0.1; %ms
% revId = fopen('HH_1D.rev', 'w');
% fprintf(revId, '<Mapping Type="Reversal">\n');
% 
% outId = fopen('HH_1D.mesh', 'w');
% fprintf(outId, 'ignore\n');
% fprintf(outId, '%f\n', timestep/1000);
% 
% formatSpec = '%.12f ';
% strip = 0;
% 
% I = 0;
% tspan = 0:timestep:10.0;
% v_0 = -52.8;
% [t,s] = ode23s(@HHv_1D, tspan, [v_0]);
% 
% gen_strip(s);
% hold on;
% 
% fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1.0);
% 
% v_0 = -52.5;
% [t,s] = ode23s(@HHv_1D, tspan, [v_0]);
% 
% gen_strip(s);
% 
% v_0 = -100;
% [t,s] = ode23s(@HHv_1D, tspan, [v_0]);
% 
% vs = s(:,1);
% 
% gen_strip(s);
% 
% fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1.0);
% 
% fprintf(outId, 'end\n');
% 
% fprintf(revId, '</Mapping>\n');
% fclose(revId);
% 
% outId = fopen('HH_1D.stat', 'w');
% fprintf(outId, '<Stationary>\n');
% fprintf(outId, '<Quadrilateral><vline>%f %f %f %f</vline><wline>%f %f %f %f</wline></Quadrilateral>\n',-64.6, -64.6, -64.2, -64.2, 0.0, 0.05, 0.05, 0.0);
% fprintf(outId, '<Quadrilateral><vline>%f %f %f %f</vline><wline>%f %f %f %f</wline></Quadrilateral>\n',-52.8, -52.8, -52.5, -52.5, 0.0, 0.05, 0.05, 0.0);
% fprintf(outId, '</Stationary>\n');
% fclose(outId);

function a = alpha_m(v)
    a = (2.5 - 0.1*v) ./ (exp(2.5 - 0.1*v) - 1);
end

function b = beta_m(v)
    b = 4*exp(-v./18);
end

function m = m_inf(v)
    a = alpha_m(v);
    b = beta_m(v);
    m = a ./ (a + b);
end

function m = tau_m(v)
    a = alpha_m(v);
    b = beta_m(v);
    m = 1 ./ (a + b);
end

function m_prime = HHm(t, vs)
    m = vs(1);
    v = vs(2);
    m_prime = - ( m - m_inf(v) ) ./ tau_m(v);
end

function a = alpha_h(v)
    a = 0.07*exp(-v ./ 20);
end

function b = beta_h(v)
    b = 1 ./ (exp(3 - 0.1*v) + 1);
end

function h = h_inf(v)
    a = alpha_h(v);
    b = beta_h(v);
    h = a ./ (a + b);
end

function h = tau_h(v)
    a = alpha_h(v);
    b = beta_h(v);
    h = 1 ./ (a + b);
end

function h_prime = HHh(t, vs)
    h = vs(1);
    v = vs(2);
    h_prime = - ( h - h_inf(v) ) ./ tau_h(v);
end

function a = alpha_n(v)
    a = (0.1 - 0.01*v) ./ (exp(1 - 0.1*v) - 1);
end

function b = beta_n(v)
    b = 0.125*exp(-v./80);
end

function n = n_inf(v)
    a = alpha_n(v);
    b = beta_n(v);
    n = a ./ (a + b);
end

function n = tau_n(v)
    a = alpha_n(v);
    b = beta_n(v);
    n = 1 ./ (a + b);
end

function n_prime = HHn(t, vs)
    n = vs(1);
    v = vs(2);
    n_prime = - ( n - n_inf(v) ) ./ tau_n(v);
end

function z_prime = HHz(t, vs)
    global alpha;
    global v_rest;
    
    z = vs(1);
    v = vs(2);
    z_prime = (-cos(alpha) * (z*cos(alpha) + n_inf(v_rest) - n_inf(v)) / tau_n(v)) - (sin(alpha) * (z*sin(alpha) + h_inf(v_rest) - h_inf(v)) / tau_h(v));
end

function w = HHw(z)
    global alpha;
    global v_rest;
    
    w = (-tan(alpha)*n_inf(v_rest)) - (z*sin(alpha));
end

function vec = HHv_full(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    
    v = vs(1);
    m = vs(2);
    h = vs(3);
    n = vs(4);
    
    m_prime = HHm(t,[m v]);
    h_prime = HHh(t,[h v]);
    n_prime = HHn(t,[n v]);
    v_prime = ( 1.0 / C ) * (  -(g_na * (m^3) * h * (v - v_na)) - (g_k * (n^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;m_prime;h_prime;n_prime];
end

function vec = HHv_m_inf(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    
    v = vs(1);
    h = vs(3);
    n = vs(4);
    
    h_prime = HHh(t,[h v]);
    n_prime = HHn(t,[n v]);
    v_prime = ( 1.0 / C ) * (  -(g_na * (m_inf(v)^3) * h * (v - v_na)) - (g_k * (n^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;h_prime;n_prime];
end

function vec = HHv_nh(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    
    v = vs(1);
    m = vs(2);
    h = vs(3);
    
    m_prime = HHm(t,[m v]);
    h_prime = HHh(t,[h v]);
    v_prime = ( 1.0 / C ) * (  -(g_na * (m^3) * h * (v - v_na)) - (g_k * ((0.79 - (1.06*h))^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;m_prime;h_prime];
end

function vec = HHv_2D_nh_m_inf(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    
    v = vs(1);
    h = vs(2);
    
    k = 0.79;
    
    h_prime = HHh(t,[h v]);
    v_prime = ( 1.0 / C ) * ( - (g_na * m_inf(v)^3 * h * (v - v_na)) - (g_k * ((k - (1.06*h))^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;h_prime];
end

function vec = HHv_1D_nh_m_inf(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    
    v = vs(1);
    
    k = 0.79;
    
    v_prime = ( 1.0 / C ) * ( - (g_na * m_inf(v)^3 * h_inf(v) * (v - v_na)) - (g_k * ((k - (1.06*h_inf(v)))^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime];
end

function [value, isterminal, direction] = HHv_reset(t, vs)
    value      = (vs(1) > 30)-0.5;
    isterminal = 1;
    direction  = 0;
end

function vec = HHv_1D(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    global v_rest;
    global alpha;
    global z_0;
    
    v = vs(1);
    
    h_0 = h_inf(-64);
    z_0 = (h_0 - h_inf(-80)) / sin(alpha);
    n = n_inf(v_rest) + z_0*cos(alpha);
    h = h_inf(v_rest) + z_0*sin(alpha);
    v_prime = ( 1.0 / C ) * ( - (g_na * m_inf(v)^3 * h * (v - v_na)) - (g_k * n^4 * (v - v_k)) - (g_l * (v - v_l)) + I);
    
    vec = [v_prime];
end

function a = gen_strip(s)

    global timestep;
    global revId;
    global outId;
    global strip;
    global formatSpec;

    right_curve = [];
    left_curve = [];

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];

    prev = s(1:1);
    for i = 2:length(s)
        p_left = [s(i,1),0.0];
        p_right = [s(i,1),0.0] + [0,0.05];
        
        if abs(prev - s(i,1)) < 0.001
            break;
        end
        
        prev = s(i,1);
        
        right_curve = [right_curve;p_right];
        left_curve = [left_curve;p_left];

        a = p_left(1);
        b = p_left(2);
        c = p_right(1);
        d = p_right(2);

        svs_1 = [svs_1, a];
        sus_1 = [sus_1, b];

        svs_2 = [svs_2, c];
        sus_2 = [sus_2, d];

        plot([a c],[b d],'k');
            hold on;
    end
    % 
    stat_a = p_right;
    stat_d = p_left;

    plot(right_curve(:,1),right_curve(:,2),'k');
    plot(left_curve(:,1),left_curve(:,2),'k');

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');

    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');

    strip = strip + 1;

    fprintf(outId, 'closed\n');
    
    a = 1;
end



