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
v_na = 55;
g_k = 100;
v_k = -80;
g_l = 0.51;
v_l = -64;

v_0 = -64;
h_0 = h_inf(-64);
m_0 = m_inf(-64);
n_0 = n_inf(-64);

% Input current
I = 4;

% v_rest is the resting potential of this neuron
% User must find the resting potential :

global approx_a;
global approx_b;
global alpha;
global v_rest;

v_rest = -64; %-45

% alpha is the angle of the slope of the curve (n_inf, h_inf) at v_rest
% User must find the angle:

alpha = -1.37; %-0.785

approx_a = -tan(alpha);
approx_b = (approx_a * n_inf(v_rest)) + h_inf(v_rest);

vs = -100:0.01:10;
hs = h_inf(vs);
ns = n_inf(vs);
plot(ns,hs,'k');
hold on;

zs = 0:0.01:1;
hs = h_inf(v_rest) + zs*sin(alpha); 
ns = n_inf(v_rest) + zs*cos(alpha); 
plot(ns,hs,'r');

ns = 0:0.01:1;
hs = approx_b - approx_a * ns;
plot(ns,hs,'b');


z_0 = (h_0 - h_inf(v_rest)) / sin(alpha)

% Plot time series of full model and 2D model
% tspan = 0:0.001:50.0;
% [t,s] = ode23s(@HHv_2D, tspan, [v_0 z_0]);
% vs = s(:,1);
% plot(t,vs,'k');
% hold on;
% 
% [t,s] = ode23s(@HHv_full, tspan, [v_0 m_0 h_0 n_0]);
% vs = s(:,1);
% plot(t,vs,'r');
% 
% hs = s(:,3);
% ns = s(:,4);
% plot(t,hs,'g');
% plot(t,ns,'b');

function m = m_inf(v)
    % to be completed by user
    % example from Rybak :
    m = 1.0 / (1 + exp(-(v + 35) / 7.8));
end

function m = tau_m(v)
    % to be completed by user
    % example from Rybak :
    m = 1;
end

function m_prime = HHm(t, vs)
    m = vs(1);
    v = vs(2);
    % to be completed by user
    % example from Rybak :
    m_prime = ( m_inf(v) - m ) / tau_m(v);
end

function h = h_inf(v)
    % to be completed by user
    % example from Rybak :
    h = 1.0 ./ (1 + exp((v + 55) ./ 7));
end

function h = tau_h(v)
    % to be completed by user
    % example from Rybak :
    h = 30 / ( (exp((v + 50)/15)) + (exp(-(v + 50)/16)) );
end

function h_prime = HHh(t, vs)
    h = vs(1);
    v = vs(2);
    % to be completed by user
    % example from Rybak :
    h_prime = ( h_inf(v) - h ) / tau_h(v);
end

function n = n_inf(v)
    % to be completed by user
    % example from Rybak :
    n = 1.0 ./ (1 + exp(-(v + 28) ./ 15));
end

function n = tau_n(v)
    % to be completed by user
    % example from Rybak :
    n = 7 / ( (exp((v + 40)/40)) + (exp(-(v + 40)/50)) );
end

function n_prime = HHn(t, vs)
    n = vs(1);
    v = vs(2);
    % to be completed by user
    % example from Rybak :
    n_prime = ( n_inf(v) - n ) / tau_n(v);
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
    v_prime = ( 1.0 / C ) * (  -(g_na * (m_inf(v)^3) * h * (v - v_na)) - (g_k * (n^4) * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;m_prime;h_prime;n_prime];
end

function vec = HHv_2D(t, vs)
    global C;
    global g_na;
    global v_na;
    global g_k;
    global v_k;
    global g_l;
    global v_l;
    global I;
    global approx_a;
    global approx_b;
    global v_rest;
    global alpha;
    
    v = vs(1);
    z = vs(2);
    
    z_prime = HHz(t,[z v]);
    w = HHw(z);
    n = n_inf(v_rest) + z*cos(alpha);
    h = h_inf(v_rest) + z*sin(alpha);
    v_prime = ( 1.0 / C ) * ( - (g_na * m_inf(v)^3 * h * (v - v_na)) - (g_k * n^4 * (v - v_k)) - (g_l * (v - v_l)) + I);
    vec = [v_prime;z_prime];
end




