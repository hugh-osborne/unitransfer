%#ok<*NOPTS>
% Reset should be -50, threshold is -40,  d = 0 (for now)
timestep = 1;
revId = fopen('rinzel.rev', 'w');
fprintf(revId, '<Mapping Type="Reversal">\n');

outId = fopen('rinzel.mesh', 'w');
fprintf(outId, 'ignore\n');
fprintf(outId, '%f\n', timestep/1000);

formatSpec = '%.12f ';
strip = 1;

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

stat_a = 0;
stat_b = 0;
stat_c = 0;
stat_d = 0;

centre_point_x = -52.85;

centre_point_y_hi = ((I - g_l.*(centre_point_x-E_l)).*(1 + exp((centre_point_x-theta_m)/sig_m)))./(g_nap.*(centre_point_x-E_na)) + 0.001;
centre_point_y_lo = ((I - g_l.*(centre_point_x-E_l)).*(1 + exp((centre_point_x-theta_m)/sig_m)))./(g_nap.*(centre_point_x-E_na)) - 0.001;

start_point_x = -64.4;

start_point_y = ((I - g_l.*(start_point_x-E_l)).*(1 + exp((start_point_x-theta_m)/sig_m)))./(g_nap.*(start_point_x-E_na));

full_time = 7000;

strip_cell_strip_width = 25;

tspan = 0:timestep:full_time;
[t1,s1] = ode23s(@izshikevich, tspan, [start_point_x start_point_y]);

x1 = s1(:,1);
y1 = s1(:,2);

time_end = full_time;

if find(x1>centre_point_x,1,'first') > 1
    time_end = find(x1>centre_point_x,1,'first')*timestep
end

% Nullcline strip (lower)
right_curve = [];
left_curve = [];

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];

nullcline_strip_length = (time_end/timestep)-100;
nullcline_strip_length = nullcline_strip_length+6;

for i = 2:nullcline_strip_length
    fwd_vec = s1(i+1,:)-s1(i,:);
    right_vec = [fwd_vec(2),-fwd_vec(1)];
    right_vec =  0.01 * right_vec / norm(right_vec);
    
    left_vec = -right_vec * 0.5;
    
    p_left = s1(i,:) + [-0.005,0.0];
    p_right = s1(i,:) + [0.005,0.0];
    
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

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1.0);

strip = strip + 1;

fprintf(outId, 'closed\n');

% (LEFT SECTION) Starting from Nullcline strip (left) going backwards to the left

v0 = left_curve(5+(1*strip_cell_strip_width),1);
h0 = left_curve(5+(1*strip_cell_strip_width),2);
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);

for i = 5+(1*strip_cell_strip_width):strip_cell_strip_width:nullcline_strip_length-(1*strip_cell_strip_width)
    v1 = left_curve(i+strip_cell_strip_width,1);
    h1 = left_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:500;

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;
    
    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end


    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 1,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, 1,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,1,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

strip_cell_strip_width = 100;

%(LOWER SECTION) Starting from the nullcline strip (right) and going
%backwards to the right.

strip_cell_strip_width = 200;

small_strip_width = strip_cell_strip_width / 20;
small_strip_start = 15;

v0 = right_curve(5+(small_strip_width),1);
h0 = right_curve(5+(small_strip_width),2);
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);

for i = 5+(small_strip_width):small_strip_width:((small_strip_start * strip_cell_strip_width)-small_strip_width)+5
    v1 = right_curve(i+small_strip_width,1);
    h1 = right_curve(i+small_strip_width,2);
    tspan = 0:timestep:150; %170 %- (30 * (max(0,(((v1-right_curve(5,1)) / (right_curve((small_strip_start * strip_cell_strip_width)-small_strip_width,1) - right_curve(5,1)))))^2));

    if strip == 531
        tspan = 0:timestep:500;
    end
    
    if strip == 352
        tspan = 0:timestep:500;
        [t0,right_last_row_1] = ode23s(@izshikevich_backward, 0:timestep:500, [v1 h1]);
    end
    
    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    
    cut = find(x1>10,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2>10,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');

    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 1,i-1, 1.0/small_strip_width);
    for j = 1:small_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, 1,i+j-1, 1.0/small_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,1,i+small_strip_width-1-1, 1.0/small_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

for i = 5+(small_strip_start*strip_cell_strip_width):strip_cell_strip_width:5+nullcline_strip_length-strip_cell_strip_width
    v1 = right_curve(i+strip_cell_strip_width,1);
    h1 = right_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:150; %160 %+ (50 * (max(0,(((v1-right_curve(((small_strip_start * strip_cell_strip_width)-small_strip_width),1)) / (right_curve((nullcline_strip_length-strip_cell_strip_width),1)-right_curve(((small_strip_start * strip_cell_strip_width)-small_strip_width),1)))))^2));

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    
    cut = find(x1>10,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2>10,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');

    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 1,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, 1,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,1,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

i = i + strip_cell_strip_width;

strip_cell_strip_width = 100;

v1 = right_curve(i+strip_cell_strip_width,1);
    h1 = right_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:150; %+ (50 * (max(0,(((v1-right_curve(((small_strip_start * strip_cell_strip_width)-small_strip_width),1)) / (right_curve((nullcline_strip_length-strip_cell_strip_width),1)-right_curve(((small_strip_start * strip_cell_strip_width)-small_strip_width),1)))))^2));

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    
    cut = find(x1>10,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2>10,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');

    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 1,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, 1,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,1,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
 
    
strip_cell_strip_width = 100;

% % Nullcline strip (upper)
% 
start_point_x = -61.255;

start_point_y = ((I - g_l.*(start_point_x-E_l)).*(1 + exp((start_point_x-theta_m)/sig_m)))./(g_nap.*(start_point_x-E_na));

full_time = 3920;

tspan = 80:timestep:full_time;
[t1,s1] = ode23s(@izshikevich, tspan, [start_point_x start_point_y]);

x1 = s1(:,1);
y1 = s1(:,2);

time_end = full_time;

right_curve = [];
left_curve = [];

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];

nullcline_strip_length = (time_end/timestep)-100;
nullcline_strip_length = nullcline_strip_length + 6;

nullcline_strip_num = strip;

for i = 2:nullcline_strip_length
    fwd_vec = s1(i+1,:)-s1(i,:);
    right_vec = [fwd_vec(2),-fwd_vec(1)];
    right_vec =  0.01 * right_vec / norm(right_vec);
    
    left_vec = -right_vec;
    
    p_left = s1(i,:) + [0.005 + ((1-(i/nullcline_strip_length))*0.05),0.0];
    p_right = s1(i,:) + [-0.005 - ((1-(i/nullcline_strip_length))*0.05),0.0];
    
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

stat_b = p_left;
stat_c = p_right;

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

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1.0);

strip = strip + 1;

fprintf(outId, 'closed\n');

% (LEFT SECTION) Starting from Nullcline strip (left) going backwards to the left
v0 = right_curve(5,1);
h0 = right_curve(5,2);
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);
small_width = strip_cell_strip_width / 10;

left_last_row_1 = prev_s2;

for i = 5:small_width:(4*strip_cell_strip_width)-small_width
    v1 = right_curve(i+small_width,1);
    h1 = right_curve(i+small_width,2);
    tspan = 0:timestep:500;

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/small_width);
    for j = 1:small_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/small_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+small_width-1-1, 1.0/small_width);
   
    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

for i = (4*strip_cell_strip_width):strip_cell_strip_width:nullcline_strip_length-(strip_cell_strip_width)
    v1 = right_curve(i+strip_cell_strip_width,1);
    h1 = right_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:500;

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);
   
    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

final_strip_width = nullcline_strip_length - (i+strip_cell_strip_width);

v1 = right_curve(nullcline_strip_length-1,1);
h1 = right_curve(nullcline_strip_length-1,2);
tspan = 0:timestep:500;

s1 = prev_s2;
[t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
prev_s2 = s2;

x1 = s1(:,1);

cut = find(x1<-120,1,'first');
if cut > 1
    x1 = x1(1:cut-1);
    y1 = s1(1:cut-1,2);
else
    x1 = s1(:,1);
    y1 = s1(:,2);
end

x2 = s2(:,1);

cut = find(x2<-120,1,'first');
if cut > 1
    x2 = x2(1:cut-1);
    y2 = s2(1:cut-1,2);
else
    x2 = s2(:,1);
    y2 = s2(:,2);
end

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 1:min(length(x2),length(x1))
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);

    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];

    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];

    plot([a c],[b d],'k');
    hold on;
end

svs_1 = fliplr(svs_1);
sus_1 = fliplr(sus_1);

svs_2 = fliplr(svs_2);
sus_2 = fliplr(sus_2);

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/final_strip_width);
for j = 1:final_strip_width-2
    fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/final_strip_width);
end
fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+final_strip_width-1-1, 1.0/final_strip_width);

strip = strip + 1;

axis fill
xlim([-70 -40]);
ylim([0 1.2]);


title('Rinzel');
ylabel('h');
xlabel('v');

plot(x1,y1,'k');
plot(x2,y2,'k');
hold on;

% (UPPER SECTION) Starting from the nullcline strip (right) and going
% backwards to the right.
reset_cell = -83+(4*strip_cell_strip_width);
v0 = left_curve(-83+(4*strip_cell_strip_width),1);
h0 = left_curve(-83+(4*strip_cell_strip_width),2);
tspan = 0:timestep:220;
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);

for i = -83+(4*strip_cell_strip_width):strip_cell_strip_width:nullcline_strip_length-strip_cell_strip_width
    v1 = left_curve(i+strip_cell_strip_width,1);
    h1 = left_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:140 + (90 * (1-(max(0,(((right_curve((5+(1*strip_cell_strip_width)),1)-v1) / (right_curve((5+(1*strip_cell_strip_width)),1)-right_curve((nullcline_strip_length-strip_cell_strip_width),1)))))^2))); %130

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;

    x1 = s1(:,1);
    y1 = s1(:,2);
    
    cut = find(x1 < -90 | y1 < 0 ,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    y2 = s2(:,2);
    
    cut = find(x2 < -90 | y2 < 0,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');

    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

final_strip_width = nullcline_strip_length - (i+strip_cell_strip_width);

v1 = left_curve(nullcline_strip_length-1,1);
h1 = left_curve(nullcline_strip_length-1,2);
tspan = 0:timestep:140; %130

s1 = prev_s2;
[t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
prev_s2 = s2;

x1 = s1(:,1);
y1 = s1(:,2);

cut = find(x1 < -90 | y1 < 0 ,1,'first');
if cut > 1
    x1 = x1(1:cut-1);
    y1 = s1(1:cut-1,2);
else
    x1 = s1(:,1);
    y1 = s1(:,2);
end

x2 = s2(:,1);
y2 = s2(:,2);

cut = find(x2 < -90 | y2 < 0,1,'first');
if cut > 1
    x2 = x2(1:cut-1);
    y2 = s2(1:cut-1,2);
else
    x2 = s2(:,1);
    y2 = s2(:,2);
end

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 1:min(length(x2),length(x1))
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);

    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];

    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];

    plot([a c],[b d],'k');
    hold on;
end

svs_1 = fliplr(svs_1);
sus_1 = fliplr(sus_1);

svs_2 = fliplr(svs_2);
sus_2 = fliplr(sus_2);

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/final_strip_width);
for j = 1:final_strip_width-2
    fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/final_strip_width);
end
fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+final_strip_width-1-1, 1.0/final_strip_width);

strip = strip + 1;

axis fill
xlim([-70 -40]);
ylim([0 1.2]);


title('Rinzel');
ylabel('h');
xlabel('v');

plot(x1,y1,'k');
plot(x2,y2,'k');
hold on;

nullcline_strip_attach = strip;

vstep = 0.25;
v_init = -61.3;
num_strips_to_join = length(v_init:0.1:-60.06);
h_init = 0.9175;

[t0,left_last_row_2] = ode45(@izshikevich_backward, 0:timestep:200, [v_init h_init]);

for v0 = v_init:0.1:-60.06
    
    cell_limit = 200 + (200 *((v0 - v_init) / (-60-v_init))) ;
    h0 = h_init;
    tspan = 0:timestep:cell_limit;
    
    [t1,s1] = ode45(@izshikevich_backward, tspan, [v0 h0]);
    [t2,s2] = ode45(@izshikevich_backward, tspan, [v0+0.1 h0]);

    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    if strip == nullcline_strip_attach
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,0, 1.0);
    else
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, strip+num_strips_to_join,0, 1.0);
    end
    
    strip = strip + 1;
    
    num_strips_to_join = num_strips_to_join + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    hold on;
end

left_last_row_3 = s2;

% Connect upper section and mid sections left_last_row_1 and left_last_row_2

% left_last_row_1 = fliplr(left_last_row_1);
% left_last_row_2 = fliplr(left_last_row_2);

x1 = left_last_row_1(:,1);
y1 = left_last_row_1(:,2);

x1 = fliplr(x1);
special_x1 = x1;
y1 = fliplr(y1);
special_y1 = y1;

x2 = left_last_row_2(:,1);
y2 = left_last_row_2(:,2);

x2 = fliplr(x2);
y2 = fliplr(y2);

special_x2 = x2;
special_y2 = y2;


svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 1:163
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);

    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];

    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];

    plot([a c],[b d],'k');
    hold on;
end

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, strip-1, length(svs_1)-2, 1.0);

strip = strip + 1;

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');


vstep = 0.25;
v_init = -61.2;
h_init = 0.9175;
for v0 = v_init:0.1:-60.06
    
    cell_limit = 280;
    h0 = h_init;
    tspan = 0:timestep:cell_limit;
    
    [t1,s1] = ode45(@izshikevich, tspan, [v0 h0]);
    [t2,s2] = ode45(@izshikevich, tspan, [v0+0.1 h0]);

    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,reset_cell, 1.0);
    
    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    hold on;
end

% Nullcline strip (right)

start_point_x = 10;

start_point_y = ((I - g_l.*(start_point_x-E_l)).*(1 + exp((start_point_x-theta_m)/sig_m)))./(g_nap.*(start_point_x-E_na)) + 0.01;

full_time = 800;

tspan = 0:timestep:full_time;
[t1,s1] = ode23s(@izshikevich, tspan, [start_point_x start_point_y]);

x1 = s1(:,1);
y1 = s1(:,2);

time_end = full_time;

right_curve = [];
left_curve = [];

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];

nullcline_strip_length = (time_end/timestep)-100;
nullcline_strip_length = nullcline_strip_length + 6;

nullcline_strip_num = strip;

for i = 1:nullcline_strip_length
    fwd_vec = s1(i+1,:)-s1(i,:);
    right_vec = [fwd_vec(2),-fwd_vec(1)];
    right_vec =  0.0005 * right_vec / norm(right_vec);
    
    left_vec = -right_vec * 0.5;
    
    p_left = s1(i,:) + left_vec + [0.0,0.0];
    p_right = s1(i,:) + right_vec + [0.0,0.0];
    
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

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 93, 11, 1.0);

strip = strip + 1;

fprintf(outId, 'closed\n');


% (RIGHT SECTION) Starting from far to the right, going backward towards the nullcline

v0 = right_curve(19,1);
h0 = right_curve(19,2);
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);

strip_cell_strip_width = 1;
for i = 20:strip_cell_strip_width:100-strip_cell_strip_width
    v1 = right_curve(i+strip_cell_strip_width,1);
    h1 = right_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:80; %100

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;
    
    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end


    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

strip_cell_strip_width = 2;
for i = 100:strip_cell_strip_width:nullcline_strip_length-strip_cell_strip_width
    v1 = right_curve(i+strip_cell_strip_width,1);
    h1 = right_curve(i+strip_cell_strip_width,2);
    tspan = 0:timestep:80 %+ (100 * (max(0,(((-20-v1) / (25))))^2)); % 70

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;
    
    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end


    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i+j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i+strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

v0 = left_curve(nullcline_strip_length-1,1);
h0 = left_curve(nullcline_strip_length-1,2);
[t1,prev_s2] = ode23s(@izshikevich_backward, tspan, [v0 h0]);

[t0,right_last_row_2] = ode23s(@izshikevich_backward, 0:timestep:200, [v0 h0]);

strip_cell_strip_width = 20;
for i = nullcline_strip_length-1:-strip_cell_strip_width:nullcline_strip_length-(30*strip_cell_strip_width)
    v1 = left_curve(i-strip_cell_strip_width,1);
    h1 = left_curve(i-strip_cell_strip_width,2);
    tspan = 0:timestep:200;

    s1 = prev_s2;
    [t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
    prev_s2 = s2;
    
    x1 = s1(:,1);
    
    cut = find(x1>10,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2>10,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end


    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end
    
    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, nullcline_strip_num,i-1, 1.0/strip_cell_strip_width);
    for j = 1:strip_cell_strip_width-2
        fprintf(revId, '%i,%i\t%i,%i\t%f\n',strip, 0, nullcline_strip_num,i-j-1, 1.0/strip_cell_strip_width);
    end
    fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0,nullcline_strip_num,i-strip_cell_strip_width-1-1, 1.0/strip_cell_strip_width);

    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    plot(x2,y2,'k');
    hold on;
end

vstep = 0.25;
v_init = 8.72;
h_init = 0.662;

strip_count = 0;
for h0 = h_init:0.01:1.2
    
    strip_count = strip_count + 1;
    
    cell_limit = 110 + (300 *((h0 - h_init) / (1.2-h_init))) ;
    
    if strip_count == 19
        cell_limit = 500;
        [t0,left_last_row_4] = ode45(@izshikevich_backward, 0:timestep:500, [v0 h0]);
    end
    
    if strip_count >= 20 && strip_count <= 21
        cell_limit = 500;
    end
    
    v0 = v_init;
    tspan = 0:timestep:cell_limit;
    
    [t1,s1] = ode45(@izshikevich_backward, tspan, [v0 h0]);
    [t2,s2] = ode45(@izshikevich_backward, tspan, [v0 h0 + 0.01]);

    x1 = s1(:,1);
    
    cut = find(x1<-120,1,'first');
    if cut > 1
        x1 = x1(1:cut-1);
        y1 = s1(1:cut-1,2);
    else
        x1 = s1(:,1);
        y1 = s1(:,2);
    end

    x2 = s2(:,1);
    
    cut = find(x2<-120,1,'first');
    if cut > 1
        x2 = x2(1:cut-1);
        y2 = s2(1:cut-1,2);
    else
        x2 = s2(:,1);
        y2 = s2(:,2);
    end

    svs_1 = [];
    sus_1 = [];
    svs_2 = [];
    sus_2 = [];
    for t = 1:min(length(x2),length(x1))
        a = x1(t);
        b = y1(t);
        c = x2(t);
        d = y2(t);

        svs_1 = [svs_1, x1(t)];
        sus_1 = [sus_1, y1(t)];
        
        svs_2 = [svs_2, x2(t)];
        sus_2 = [sus_2, y2(t)];

        plot([a c],[b d],'k');
        hold on;
    end

    svs_1 = fliplr(svs_1);
    sus_1 = fliplr(sus_1);
    
    svs_2 = fliplr(svs_2);
    sus_2 = fliplr(sus_2);

    fprintf(outId, formatSpec, svs_1);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_1);
    fprintf(outId, '\n');
    
    fprintf(outId, formatSpec, svs_2);
    fprintf(outId, '\n');
    fprintf(outId, formatSpec, sus_2);
    fprintf(outId, '\n');
    fprintf(outId, 'closed\n');
    
    strip = strip + 1;
    
    axis fill
    xlim([-70 -40]);
    ylim([0 1.2]);
    

    title('Rinzel');
    ylabel('h');
    xlabel('v');

    plot(x1,y1,'k');
    hold on;
end

% Connect upper section and mid sections left_last_row_3 and left_last_row_4

x1 = left_last_row_3(:,1);
y1 = left_last_row_3(:,2);

special_x1 = x1;
special_y1 = y1;

x2 = left_last_row_4(:,1);
y2 = left_last_row_4(:,2);


special_x2 = x2;
special_y2 = y2;


svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
diff = 293 - 275;
for t = 293:-1:73
    a = x1(t-diff);
    b = y1(t-diff);
    c = x2(t);
    d = y2(t);

    svs_1 = [svs_1, x1(t-diff)];
    sus_1 = [sus_1, y1(t-diff)];

    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];

    plot([a c],[b d],'k');
    hold on;
end

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, strip-1, length(svs_1)-2, 1.0);

strip = strip + 1;

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');

% Connect upper section and mid sections right_last_row_1 and
% right_last_row_2

x1 = right_last_row_1(:,1);
y1 = right_last_row_1(:,2);

x1 = flipud(x1);
y1 = flipud(y1);

special_x1 = x1;
special_y1 = y1;

x2 = right_last_row_2(:,1);
y2 = right_last_row_2(:,2);

x2 = flipud(x2);
y2 = flipud(y2);

special_x2 = x2;
special_y2 = y2;

svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 12:59
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);

    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];

    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];

    plot([a c],[b d],'k');
    hold on;
end

fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, strip-1, length(svs_1)-2, 1.0);

strip = strip + 1;

fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');

fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');

v0 = stat_a(1);
h0 = stat_a(2);
v1 = stat_b(1);
h1 = stat_b(2);
tspan = 0:timestep:150;
 
[t1,s1] = ode23s(@izshikevich_backward, tspan, [v0 h0]);
[t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
 
x1 = s1(:,1);
 
cut = find(x1<-120,1,'first');
if cut > 1
    x1 = x1(1:cut-1);
    y1 = s1(1:cut-1,2);
else
    x1 = s1(:,1);
    y1 = s1(:,2);
end
 
x2 = s2(:,1);
 
cut = find(x2<-120,1,'first');
if cut > 1
    x2 = x2(1:cut-1);
    y2 = s2(1:cut-1,2);
else
    x2 = s2(:,1);
    y2 = s2(:,2);
end
 
svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 1:min(length(x2),length(x1))
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);
 
    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];
 
    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];
 
    plot([a c],[b d],'k');
    hold on;
end
 
svs_1 = fliplr(svs_1);
sus_1 = fliplr(sus_1);
 
svs_2 = fliplr(svs_2);
sus_2 = fliplr(sus_2);
 
fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');
 
fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');
 
fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1);
 
strip = strip + 1;
 
axis fill
xlim([-70 -40]);
ylim([0 1.2]);
 
title('Rinzel');
ylabel('h');
xlabel('v');
 
plot(x1,y1,'k');
plot(x2,y2,'k');
hold on;
 
v0 = stat_c(1);
h0 = stat_c(2);
v1 = stat_d(1);
h1 = stat_d(2);
tspan = 0:timestep:150;
 
[t1,s1] = ode23s(@izshikevich_backward, tspan, [v0 h0]);
[t2,s2] = ode23s(@izshikevich_backward, tspan, [v1 h1]);
 
x1 = s1(:,1);
 
cut = find(x1<-120,1,'first');
if cut > 1
    x1 = x1(1:cut-1);
    y1 = s1(1:cut-1,2);
else
    x1 = s1(:,1);
    y1 = s1(:,2);
end
 
x2 = s2(:,1);
 
cut = find(x2<-120,1,'first');
if cut > 1
    x2 = x2(1:cut-1);
    y2 = s2(1:cut-1,2);
else
    x2 = s2(:,1);
    y2 = s2(:,2);
end
 
svs_1 = [];
sus_1 = [];
svs_2 = [];
sus_2 = [];
for t = 1:min(length(x2),length(x1))
    a = x1(t);
    b = y1(t);
    c = x2(t);
    d = y2(t);
 
    svs_1 = [svs_1, x1(t)];
    sus_1 = [sus_1, y1(t)];
 
    svs_2 = [svs_2, x2(t)];
    sus_2 = [sus_2, y2(t)];
 
    plot([a c],[b d],'k');
    hold on;
end
 
svs_1 = fliplr(svs_1);
sus_1 = fliplr(sus_1);
 
svs_2 = fliplr(svs_2);
sus_2 = fliplr(sus_2);
 
fprintf(outId, formatSpec, svs_1);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_1);
fprintf(outId, '\n');
 
fprintf(outId, formatSpec, svs_2);
fprintf(outId, '\n');
fprintf(outId, formatSpec, sus_2);
fprintf(outId, '\n');
fprintf(outId, 'closed\n');
 
fprintf(revId, '%i,%i\t%i,%i\t%f\n', strip, 0, 0,0, 1);
strip = strip + 1;
 
axis fill
xlim([-70 -40]);
ylim([0 1.2]);
 
title('Rinzel');
ylabel('h');
xlabel('v');
 
plot(x1,y1,'k');
plot(x2,y2,'k');
hold on;
 
fprintf(outId, 'end\n');

%nullclines
nvs = -120:0.01:0;
    
v_nulls = ((I - g_l.*(nvs-E_l)).*(1 + exp((nvs-theta_m)/sig_m)))./(g_nap.*(nvs-E_na));
h_nulls = 1 ./ (1+exp((nvs-theta_h)./sig_h));

plot(nvs,v_nulls,'r');
plot(nvs,h_nulls,'b');

%stationary quad

stat_min_v = -59.05;
stat_max_v = -58.95;
stat_min_w = 0.56;
stat_max_w = 0.54;

plot([stat_a(1) stat_b(1)],[stat_a(2) stat_b(2)],'g');
plot([stat_b(1) stat_c(1)],[stat_b(2) stat_c(2)],'g');
plot([stat_c(1) stat_d(1)],[stat_c(2) stat_d(2)],'g');
plot([stat_d(1) stat_a(1)],[stat_d(2) stat_a(2)],'g');

fprintf(revId, '</Mapping>\n');
fclose(revId);

outId = fopen('rinzel.stat', 'w');
fprintf(outId, '<Stationary>\n');
fprintf(outId, '<Quadrilateral><vline>%f %f %f %f</vline><wline>%f %f %f %f</wline></Quadrilateral>\n', stat_a(1), stat_b(1), stat_c(1), stat_d(1), stat_a(2), stat_b(2), stat_c(2), stat_d(2));
fprintf(outId, '</Stationary>\n');
fclose(outId);
