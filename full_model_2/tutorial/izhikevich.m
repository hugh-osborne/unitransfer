function[]=izhikevich(varargin)
if nargin > 0 
	eval(varargin{1}); %execute the argument 
else
    main;
end;

%--------------------------------
function[]=main
global pauseflag drawnull thresh quitflag 
global VU1 VU2 L 
tau = 0.01;
L = 5000/tau;
drawnull=1;
quitflag=0;

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
g_l = 0.01; %mS
E_l = -64.0; %mV
I = 0.0; %

figNumber = figure(1);
clf;
set(figNumber,'NumberTitle','off','doublebuffer','on',...
        'Name','Simple Model by Izhikevich (2003)',...
        'Units','normalized','toolbar','figure',...
        'Position',[0.05 0.1 0.9 0.8]);
h1=subplot(4,2,1);
set(h1,'Position',[0.05 0.75 0.27 0.2])
vtrace1=line('color','k','LineStyle','-','erase','background','xdata',[],'ydata',[],'zdata',[]);
axis([0 5000 -100 30])
title('membrane potential, v')
xlabel('time (ms)');

h4=subplot(4,2,4);
set(h4,'Position',[0.05 0.5 0.27 0.15])
vtrace2=line('color','k','LineStyle','-','erase','background','xdata',[],'ydata',[],'zdata',[]);
axis([0 5000 -100 30])
title('membrane potential, v')
xlabel('time (ms)');

% h2=subplot(4,2,3);
% set(h2,'Position',[0.05 0.5 0.27 0.15])
% utrace=line('color','k','LineStyle','-','erase','background','xdata',[],'ydata',[],'zdata',[]);
% axis([0 100 -0.2 1.0])
% title('recovery variable, u')
% xlabel('time (ms)')

h3=subplot(2,2,2);
set(h3,'Position',[0.4 0.5 0.55 0.45])
vnull=line('color','b','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
good_vnull=line('color','r','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
h_nap_null=line('color','r','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
h_na_null=line('color','b','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
m_null=line('color','b','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
axis([-100 0 -0.2 1.0]);
title('phase portrait')
xlabel('v');ylabel('u');

RG_E = burster(-64.0, h3, 17, 0.583); %burster_min = 25
RG_E.setI(0.0);
RG_F = burster(-64.0, h3, 17, 0.583); %burster_min = 25
RG_F.setI(0.0);

RG_E.attachNext(RG_F);
RG_F.attachNext(RG_E);

VU1 = [0;0]*ones(1,L);
VU2 = [0;0]*ones(1,L);
t=0;

% The main loop
while quitflag==0
   if drawnull>0
        drawnull=0;
        vv=-100:0.1:0;
      %  set(vnull,'xdata',vv,'ydata',((I - g_l.*(vv-E_l)).*(1 + exp((vv-theta_m)/sig_m)))./(g_nap.*(vv-E_na)));
        g_na = 30;
        theta_m_na = -35;
        sig_m_na = -7.8;
        E_k= -80;
        g_k = 1;
        set(good_vnull,'xdata',vv,'ydata',((g_l.*(vv-E_l)) + (g_k.*(((1./(1 + exp((vv+28)/-15)))).^4).*(vv-E_k)) + (g_na.*0.7243.*(((1./(1 + exp((vv-theta_m_na)/sig_m_na)))).^3).*(vv-E_na)) - I ) ./ ((vv-E_na).*-g_nap.*(1./(1 + exp((vv-theta_m)/sig_m)))));
        set(h_nap_null,'xdata',vv,'ydata',1 ./ (1+exp((vv-theta_h)./sig_h)));
       % set(h_na_null,'xdata',vv,'ydata',1 ./ (1+exp((vv+55)./7)));
      %  set(m_null,'xdata',vv,'ydata',1 ./ (1+exp((-vv-28)./15)));
        drawnow;
   end
   
   RG_E.update(tau, t);
   RG_F.update(tau, t);

   t=t+tau; 
   if t > 1000
       RG_F.setI(5); %burster_min = 4200
   end
   
   if t > 0
       RG_E.setI(5); %burster_min = 4200
   end
   
   RG_E.draw();
   RG_F.draw();
   
   VU1 = [[RG_E.v;RG_E.h_nap],VU1(:,1:L-1)];
   VU2 = [[RG_F.v;RG_F.h_nap],VU2(:,1:L-1)];
   
   set(vtrace1,'xdata',tau*(L:-1:1),'ydata',VU1(1,:));
   set(vtrace2,'xdata',tau*(L:-1:1),'ydata',VU2(1,:));
   drawnow;

end



