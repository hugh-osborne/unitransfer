function izsh = izshikevich_backward( t, ws )
    % TODO : rename this function.
    
    g_nap = 0.25; %mS - %250.0; %nS
    theta_m = -47.1; %mV
    sig_m = -3.1; %mV
    theta_h = -59; %mV
    sig_h = 8; %mV
    tau_h = 1200; %ms 
    E_na = 55; %mV
    C = 1; %uF - %1000; %pF
    g_l = 0.1; %mS - %100.0; %nS
    E_l = -64.0; %mV
    I = 0.0; %
    I_h = 0; %
    
    v = ws(1);
    h = ws(2);
    
    v_prime = (((-g_nap * h * (v - E_na) * ((1 + exp((v-theta_m)/sig_m))^-1)) - (g_l*(v - E_l))) / C)+I;
    h_prime = (((1 + (exp((v - theta_h)/sig_h)))^(-1)) - h ) / (tau_h/cosh((v - theta_h)/(2*sig_h))) + I_h;
    
    izsh = [-v_prime; -h_prime];
end

