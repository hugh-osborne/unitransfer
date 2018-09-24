classdef burster < handle
    %BURSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        g_nap %mS - %250.0; %nS
        g_na
        g_k
        g_syn
        theta_m %mV
        sig_m %mV
        theta_m_na %mV
        sig_m_na %mV
        theta_h %mV
        sig_h %mV
        tau_h %ms 
        E_na %mV
        E_k %mV
        C %uF - %1000; %pF
        g_l %mS
        E_l %mV
        E_synE
        E_synI
        I
        v
        h_nap
        h_na
        h_synI
        h_synE
        mlr
        m_k
        L
        VU
        VVU
        plot
        head
        next
        prev
        syn_eff
        t_ref
        t_last
        refracting
    end
    
    methods
        function obj = burster(E_l, plt, eff, h)
            obj.g_nap = 0.25; %mS - %250.0; %nS
            obj.g_na = 30;
            obj.g_k = 6; %1
            obj.g_l = 0.1; %mS
            obj.g_syn = 0.05;
            obj.E_k = -80;
            obj.E_na = 55; %mV
            obj.E_l = -64.0; %mV
            obj.mlr = 0;  
            obj.C = 1; %uF - %1000; %pF
            % Depending on how much noise we expect, E_l should be reduced to the point
            % where there is still no activity - that is, the stationary point comes before the v
            % nullcline peak. (higher I => lower v nullcline peak - stationary point moves to the right, 
            % lower E_l => higher nullcline peak - stationary point moves to the left)
            % ***Note for each -1 in E_l, I must increase by g_l to keep the same stationary point.*** 
            obj.I = 0.0; %
            obj.v=-70;
            obj.h_nap=h;
            obj.h_na=0;
            obj.m_k=0;
            obj.L = 500;
            obj.VU = [obj.v;obj.h_nap]*ones(1,obj.L);
            obj.plot = plt;
            get(obj.plot,'Position'); 
            obj.head = line('color','r','Marker','.','markersize',20,'erase','xor','xdata',[],'ydata',[],'zdata',[]);
            obj.next = [];
            obj.syn_eff = eff;
            obj.t_last = 0;
            obj.t_ref = 5; %18
            obj.refracting = 0;
        end
        function r = attachNext(obj, n)
            obj.next = [obj.next, n];
            r = 0;
        end
        function r = setI(obj, i)
            obj.mlr = i;
            r = 0;
        end
        function r = update(obj,dt, t)
            I_nap = -obj.g_nap * obj.h_nap * (obj.v - obj.E_na) * ((1 + exp(-(obj.v+47.1)/3.1))^-1);
            I_l = -obj.g_l*(obj.v - obj.E_l);
            I_na = -obj.g_na * obj.h_na * (obj.v - obj.E_na) * (((1 + exp(-(obj.v+35)/7.8))^-1)^3);
            I_k = -obj.g_k * ((obj.m_k)^4) * (obj.v - obj.E_k);
            
%             if obj.refracting == 0
            obj.v=obj.v+dt*(I_nap + I_na + I_k + I_l + obj.mlr);
%             end

            part_1 = ((1 + (exp((obj.v + 59)/8)))^(-1)) - obj.h_nap;
            part_2 = (1200/cosh((obj.v + 59)/(16)));
            obj.h_nap=obj.h_nap+dt*(part_1 / part_2);
            
            part_1 = ((1 + (exp((obj.v + 55)/7)))^(-1)) - obj.h_na;
            part_2 = 30 / (exp((obj.v+50)/15) + exp(-(obj.v+50)/16));
            obj.h_na=obj.h_na+dt*(part_1 / part_2);
            
            part_1 = ((1 + (exp(-(obj.v + 28)/15)))^(-1)) - obj.m_k;
            part_2 = 7 / (exp((obj.v+40)/40) + exp(-(obj.v + 40)/50));
            obj.m_k=obj.m_k+dt*(part_1 / part_2);
            
%             if obj.v>-10 && obj.refracting == 0
%                 obj.t_last = t;
%                 obj.refracting = 1; 
% %                 obj.v=-50;
% %                 obj.h_nap=obj.h_nap-0.004;
%                 for n = obj.next
% %                     n.v = n.v - obj.syn_eff;
%                     n.refracting = 0;
%                 end
%             end    
%             if obj.refracting
%                 if t - obj.t_last > obj.t_ref
%                     obj.refracting = 0;
%                 end
%             end
            obj.VU = [[obj.v;obj.h_nap],obj.VU(:,1:obj.L-1)];
            r = 0;
        end
        function r = draw(obj)
            set(obj.head,'xdata',obj.v,'ydata',obj.h_nap);
            drawnow;
            r = 0;
        end
    end
    
end

