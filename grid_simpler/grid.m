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
g_na = 30;
theta_m_na = -35;
sig_m_na = -7.8;
E_k= -80;
g_k = 1;

timestep = 1;
transform = 1;
filename = 'grid_transform';

revId = fopen(strcat(filename,'.rev'), 'w');
fprintf(revId, '<Mapping Type="Reversal">\n');

outId = fopen(strcat(filename,'.mesh'), 'w');
fprintf(outId, 'ignore\n');
fprintf(outId, '%f\n', timestep / 1000);

formatSpec = '%.12f ';

axis fill
xlim([-90 -10]);
ylim([-0.4 1.2]);

title('Grid');
ylabel('h');
xlabel('v');

columns = 100;
rows = 80;

for i = 0:(1/columns):1
    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    
    for j = 0:(1/rows):1     
        x1 = (i * 50) - 90;
        y1 = (j * 1.4) - 0.4;
        
        if transform == 0
            svs_1 = [svs_1, x1];
            sus_1 = [sus_1, y1];

            svs_2 = [svs_2, x1+(0.01*50)];
            sus_2 = [sus_2, y1];
            
            plot([x1 x1+((1/columns)*50)],[y1 y1],'k');
            hold on;
        else
            t_r1 = rybak([x1,y1],timestep);
            t_r2 = rybak([x1+((1/columns)*50),y1],timestep);
            t_x1 = t_r1(1);
            t_y1 = t_r1(2);
            t_x2 = t_r2(1);
            t_y2 = t_r2(2);

            svs_1 = [svs_1, t_x1];
            sus_1 = [sus_1, t_y1];

            svs_2 = [svs_2, t_x2];
            sus_2 = [sus_2, t_y2];
            
            plot([t_x1 t_x2],[t_y1 t_y2],'k');
            hold on;
        end

    end
    
    plot(svs_1,sus_1,'k');
    hold on;

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');

    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');

    fprintf(outId, 'closed\n');

end

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];

for j = 0:(1/rows):1     
    x1 = -40;
    y1 = (j * 1.4) - 0.4;

    if transform == 0
        svs_1 = [svs_1, x1];
        sus_1 = [sus_1, y1];

        svs_2 = [svs_2, x1+40];
        sus_2 = [sus_2, y1];

        plot([x1 x1+40],[y1 y1],'k');
        hold on;
    else
        t_r1 = rybak([x1,y1],timestep);
        t_r2 = rybak([x1+40,y1],timestep);
        t_x1 = t_r1(1);
        t_y1 = t_r1(2);
        t_x2 = t_r2(1);
        t_y2 = t_r2(2);

        svs_1 = [svs_1, t_x1];
        sus_1 = [sus_1, t_y1];

        svs_2 = [svs_2, t_x2];
        sus_2 = [sus_2, t_y2];

        plot([t_x1 t_x2],[t_y1 t_y2],'k');
        hold on;
    end
end

plot(svs_1,sus_1,'k');
hold on;

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');

fprintf(outId, 'closed\n');

fprintf(outId, 'end\n');

axis fill
xlim([-90 -40]);
ylim([-0.4 1.4]);

title('Grid');
ylabel('h');
xlabel('v');

%nullclines
vv = -120:0.01:0;

v_nulls = ((g_l.*(vv-E_l)) + (g_k.*(((1./(1 + exp((vv+28)/-15)))).^4).*(vv-E_k)) + (g_na.*0.7243.*(((1./(1 + exp((vv-theta_m_na)/sig_m_na)))).^3).*(vv-E_na)) - I ) ./ ((vv-E_na).*-g_nap.*(1./(1 + exp((vv-theta_m)/sig_m))));
h_nulls = 1 ./ (1+exp((vv-theta_h)./sig_h));

plot(vv,v_nulls,'r');
plot(vv,h_nulls,'b');

fprintf(revId, '</Mapping>\n');
fclose(revId);
fclose(outId);

%stationary

stat_a = [-61.685 0.5835];
stat_c = [-61.695 0.5835];
stat_b = [-61.685 0.583];
stat_d = [-61.695 0.583];

plot([stat_a(1) stat_b(1)],[stat_a(2) stat_b(2)],'g');
plot([stat_b(1) stat_c(1)],[stat_b(2) stat_c(2)],'g');
plot([stat_c(1) stat_d(1)],[stat_c(2) stat_d(2)],'g');
plot([stat_d(1) stat_a(1)],[stat_d(2) stat_a(2)],'g');

outId = fopen(strcat(filename,'.stat'), 'w');
fprintf(outId, '<Stationary>\n');
fprintf(outId, '<Quadrilateral><vline>%f %f %f %f</vline><wline>%f %f %f %f</wline></Quadrilateral>\n', stat_a(1), stat_b(1), stat_c(1), stat_d(1), stat_a(2), stat_b(2), stat_c(2), stat_d(2));
fprintf(outId, '</Stationary>\n');
fclose(outId);

function r = rybak(ws, dt)
    % TODO : rename this function.
    
    g_nap = 0.25; %mS - %250.0; %nS
    g_na = 30;
    g_k = 1;
    theta_m = -47.1; %mV
    sig_m = -3.1; %mV
    theta_h = -59; %mV
    sig_h = 8; %mV
    tau_h = 1200; %ms 
    E_na = 55; %mV
    E_k = -80;
    C = 1; %uF - %1000; %pF
    g_l = 0.1; %mS - %100.0; %nS
    E_l = -64.0; %mV
    I = 0.0; %
    I_h = 0; %
    
    v = ws(1);
    h = ws(2);
    
    I_nap = -g_nap * h * (v - E_na) * ((1 + exp((-v-47.1)/3.1))^-1);
    I_l = -g_l*(v - E_l);
    I_na = -g_na * 0.7243 * (v - E_na) * (((1 + exp((v+35)/-7.8))^-1)^3);
    I_k = -g_k * (v - E_k) * (((1 + exp((v+28)/-15))^-1)^4);
    
    v_prime = ((I_nap + I_l + I_na + I_k) / C)+I;
    h_prime = (((1 + (exp((v - theta_h)/sig_h)))^(-1)) - h ) / (tau_h/cosh((v - theta_h)/(2*sig_h))) + I_h;
    
    r = [v + (dt*v_prime), h + dt*(h_prime)];
end

function r = rybak_backward(ws, dt)
    % TODO : rename this function.
    
    g_nap = 0.25; %mS - %250.0; %nS
    g_na = 30;
    g_k = 1;
    theta_m = -47.1; %mV
    sig_m = -3.1; %mV
    theta_h = -59; %mV
    sig_h = 8; %mV
    tau_h = 1200; %ms 
    E_na = 55; %mV
    E_k = -80;
    C = 1; %uF - %1000; %pF
    g_l = 0.1; %mS - %100.0; %nS
    E_l = -64.0; %mV
    I = 0.0; %
    I_h = 0; %
    
    v = ws(1);
    h = ws(2);
    
    I_nap = -g_nap * h * (v - E_na) * ((1 + exp((-v-47.1)/3.1))^-1);
    I_l = -g_l*(v - E_l);
    I_na = -g_na * 0.7243 * (v - E_na) * (((1 + exp((v+35)/-7.8))^-1)^3);
    I_k = -g_k * (v - E_k) * (((1 + exp((v+28)/-15))^-1)^4);
    
    v_prime = ((I_nap + I_l + I_na + I_k) / C)+I;
    h_prime = (((1 + (exp((v - theta_h)/sig_h)))^(-1)) - h ) / (tau_h/cosh((v - theta_h)/(2*sig_h))) + I_h;
    
    r = [v - (dt*v_prime), h - dt*(h_prime)];
end
