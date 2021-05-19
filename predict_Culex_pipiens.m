function [prdData, info] = predict_Culex_pipiens(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  % compute temperature correction factors
  TC = tempcorr(temp.ab, T_ref, T_A);     % most data at 28 C
  TC_WJO = tempcorr(temp.WJO, T_ref, T_A);
  
  % zero-variate data

  % life cycle
  info = 1;
  % cannot use get_tj_hex because f is not constant
  %pars_tj = [g k v_Hb v_He s_j kap kap_V];
  %[t_j, t_e, t_b, l_j, l_e, l_b, rho_j, v_Rj, u_Ee, info] = get_tj_hex(pars_tj, f);
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  E_0 = 1e3 * p_Am * initial_scaled_reserve(f, pars_UE0); % mJ, initial energy content
  Ww_0 = w_E * E_0/ mu_E/ d_E;                            % mg, initial wet weight
    
  % birth
  [t_b, l_b] = get_tb([g k v_Hb], f);
  aT_b = t_0 + t_b/ k_M/ TC;          % d, age at birth at f and T
  L_b = L_m * l_b;                    % cm, structural length at birth
  Lw_b = 10 * L_b/ del_M;             % mm, head length at birth
  Ww_b = 1e3 * L_b^3 * (1 + f * w);   % mg, wet weight at birth

  % metam
  p_J = k_J * E_Hb; pT_J = TC * p_J;% J/d, maturity maint rate during larval stages
  L_3 = L_b * prod([s_1 s_2 s_3].^0.5); % cm, structural length at instar 4
  r = k_M * g * (f/ l_b - 1)/ (f + g); % 1/d, specific growth rate
  E_R_ref = (1 - kap) * E_m * (g + l_b)/ (1 - l_b); % J/cm^3, reference value for [E_R]
  E_Rj = s_j * E_R_ref; % J/cm^3, threshold reprod buffer density at pupation
  aELE_0 = [0; f * E_m * L_b^3; L_b; 0]; % initial state vector
  [ER aELE] = ode45(@dget_aELE, [0; E_Rj], aELE_0, [], f, k_M, v/ L_b, p_J, p_Am/ L_b, E_m, g, kap, L_3, f_4);
  t_j = aELE(end,1); tT_j = t_j/ TC;    % d, age at pupation
  E = aELE(end,2);                  % J, reserve at pupation
  L_j = aELE(end,3);                % cm, structural length at pupation
  E_R = aELE(end,4);                % J, reprod buffer at pupation
  E_j = E + E_R;                    % J, reserve at pupation
  Ww_j = 1e3 * (L_j^3 + E_j * w_E/ mu_E/ d_E); % mg, wet weight at pupation
  Lw_j = 10 * L_j/ del_M;           % mm, length at metam
  s_M = L_j/ L_b;                   % -, acceleration factor
  
  % emergence
  v_j = v * s_M;         % cm/d, energy conductance of imago
  L_mj = L_m * s_M;      % cm, max length of mago
  [H aVEL] = ode45(@dget_aVEL, [0; E_He], [0; L_j^3; E + E_R; 1e-4], [], k_J, v/ L_b, v_j, E_m, L_mj, kap, g, E_G, kap_V);
  t_e = aVEL(end,1);     % d, time since pupation at emergence
  V_e = aVEL(end,2);     % cm^3, larval structure at emergence (should be close to zero)
  E_e = aVEL(end,3);     % J, reserve at emergence
  L_e = aVEL(end,4);     % cm, length at emergence
  tT_e = t_e/ TC;        % d, age at emergence
  E_e = E_e + kap_V * E_G * V_e; % J, convert last bit of larval structure to reserve instantaneously
  Ww_e = 1e3 * (L_e^3 + E_e * w_E/ mu_E/ d_E);   % mg, wet weight at emergence
  
  % instar age, length
  Lw_1 = Lw_b; Lw_2 = Lw_1 * s_1^0.5; Lw_3 = Lw_2 * s_2^0.5; Lw_4 = Lw_3 * s_3^0.5; % mm, head length of instar 1-4
  rT_j = TC * g * k_M * (f/ l_b - 1)/ (f + g); rT_4 = TC * g * k_M * (f_4/ l_b - 1)/ (f_4 + g); % 1/d, specific growth rate
  tT_1 = log(Lw_2/ Lw_1) * 3/ rT_j;  tT_2 = log(Lw_3/ Lw_2) * 3/ rT_j;
  tT_3 = log(Lw_4/ Lw_3) * 3/ rT_j;  tT_4 = log(Lw_j/ Lw_4) * 3/ rT_4;

  % reproduction (no kappa-rule)
  RT_i = 1e3 * kap_R * TC * (p_Am * s_M * L_e^2 - p_M * L_e^3 - k_J * E_He)/ E_0; % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;                  % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tj = tT_j;
  prdData.te = tT_e;
  prdData.am = aT_m;
  prdData.t1 = tT_1;
  prdData.t2 = tT_2;
  prdData.t3 = tT_3;
  prdData.t4 = tT_4;
  prdData.L1 = Lw_1;
  prdData.L2 = Lw_2;
  prdData.L3 = Lw_3;
  prdData.L4 = Lw_4;
  prdData.Ww0 = Ww_0;
  prdData.Wwj = Ww_j;
  prdData.Wwe = Ww_e;
  prdData.E0 = E_0;
  prdData.Ri = RT_i;
  
  % uni-variate data
  
  % time-weight 
  L_3 = L_b * prod([s_1 s_2 s_3].^0.5); % cm, structural length at instar 4
  ELE_0 = [f * E_m * L_b^3; L_b; 0];    % initial state vector
  [t ELE] = ode45(@dget_ELE, tW(:,1), ELE_0, [], f_tW, TC * k_M, TC * v/ L_b, TC * p_J, TC * p_Am/ L_b, E_m, g, kap, L_3, f_4);
  EWw = 1e3 * (ELE(:,2).^3 + sum(ELE(:,[1 3]), 2)/ mu_E * w_E/ d_E); % mg, wet weight
  
  % WJO-data
  % larva
  [W ELE] = ode45(@dget_ELE_W, [Ww_b; WJO(1:3,1)], ELE_0, [], f, k_M, v/ L_b, p_J, p_Am/L_b, E_m, g, kap, L_3, f_4, w_E, mu_E, d_V);
  ELE(1,:) = []; L = ELE(:,2); V = L.^3;  E = ELE(:,1); e = E ./ V/ E_m;  % -, scaled reserve density
  r = g * k_M * (e/ l_b - 1) ./ (e + g); % 1/d, specific growth rate
  pA = f * p_Am/ L_b * V; pD = p_J + p_M * V; %pA = 0 * pA;
  p_C = E .* (v/ L_b - r); pG = kap * p_C - p_M * V;  
  pADG = [pA, pD, pG]; % scaled growth (page 78 DEB book)
  JM = - (n_M\n_O) * eta_O * pADG'; % mol/d, mineral fluxes
  JO = - JM(3,:)'/ 24;       % mol O2/h, O2-consumption 
  EJO_l = TC_WJO(1:3) .* JO * 24.4e6 ./ WJO(1:3,1); % mm^3/g/h, specific O2 consumption 
  % pupa
  v_j = v * s_M;    % cm/d, energy conductance of imago
  L_mj = L_m * L_j/ L_b; % cm, max length of imago
  HVEL_0 = [0; L_j^3; E_j; 1e-4];
  [W HVEL] = ode45(@dget_HVEL_W, [Wwj; WJO(4:5,1)], HVEL_0, [], k_J, v/ L_b, v_j, E_m, L_mj, kap, g, mu_E, kap_V, w_E, d_E);
  HVEL(1,:) = [];
  H = HVEL(:,1); V_l = HVEL(:,2); E = HVEL(:,3); L = HVEL(:,4); V = L.^3;  e = E ./ V/ E_m;  % -, scaled reserve density
  r = v_j * (e ./ L - 1/ L_m) ./ (e + g); % 1/d, specific growth rate
  pA = [0; 0]; p_C = E .* (v_j ./ L - r); % J/d, mobilisation rate
  pD = (1 - kap) * p_C + p_M * V;  % J/d, dissipation rate
  pG = kap * p_C - p_M * V;  % J/d, change in structure
  pADG = [pA, pD, pG]; % scaled growth (page 78 DEB book)
  JM = - (n_M\n_O) * eta_O * pADG'; % mol/d, mineral fluxes
  JO = - JM(3,:)'/ 24;       % mol O2/h, O2-consumption 
  EJO_p = TC_WJO(4:5) .* JO * 24.4e6 ./ WJO(4:5,1); % mm^3/g/h, specific O2 consumption 
  % imago
  E = (WJO(6:7,1) - L_e^3) * mu_E/ w_E; % J, reserve
  e = E ./ L_e^3/ E_m;  % -, scaled reserve density
  pA = f * p_Am * s_M * L_e^2 * ones(2,1); pA = 0 * pA;
  p_C = E .* v_j/ L_e ;          % J/d, mobilisation rate
  pD = k_J * E_He + p_M * L_e^3 + ((1 - kap) * p_C - k_J * E_He) * (1 - kap_R);  % J/d, dissipation rate
  pG = 0 * pA;  % J/d, change in structure
  pADG = [pA, pD, pG]; % scaled growth (page 78 DEB book)
  JM = - (n_M\n_O) * eta_O * pADG'; % mol/d, mineral fluxes
  JO = - JM(3,:)'/ 24;       % mol O2/h, O2-consumption 
  EJO_e = TC_WJO(6:7) .* JO * 24.4e6 ./ WJO(6:7,1); % mm^3/g/h, specific O2 consumption 
  EJT_O = [EJO_l; EJO_p; EJO_e];

  % pack to output
  prdData.tW = EWw;
  prdData.WJO = EJT_O;
  
end

%% subfunctions

    function dELE = dget_ELE(t, ELE, f, k_M, k_E, p_J, p_Am, E_m, g, kap, L_3, f_4)
      % larval development
      % t: time since birth
      
      % unpack variables
      E  = ELE(1); % J, reserve
      L  = ELE(2); % cm, structural length
      %E_R= ELE(3); % J, reproduction buffer
      
      if L > L_3 % less feeding in instar 4
        f = f_4;
      end
      
      V = L^3;                          % cm^3, structural volume
      e = E/ V/ E_m;                    % -, scaled reserve density
      r = (e * k_E - g * k_M)/ (e + g); % 1/d, specific growth rate
      p_C = E * (k_E - r);              % J/d, mobilisation rate
      
      dE = f * p_Am * V - p_C;          % J/d, change in reserve
      dL = r * L/ 3;                    % cm/d, change in length
      dE_R = (1 - kap) * p_C - p_J;     % J/d, change in reprod buffer

      dELE = [dE; dL; dE_R];   
    end

    function daELE = dget_aELE(ER, aELE, f, k_M, k_E, p_J, p_Am, E_m, g, kap, L_3, f_4)
      % larval development
      % ER condition index: [E_R] = E_R/V in J/cm^3
      
      % unpack variables
      a  = aELE(1); % d, time since birth
      E  = aELE(2); % J, reserve
      L  = aELE(3); % cm, structural length
      E_R= aELE(4); % J, reproduction buffer
      
      if L > L_3 % less feeding in instar 4
        f = f_4;
      end
      
      V = L^3;                          % cm^3, structural volume
      e = E/ V/ E_m;                    % -, scaled reserve density
      r = (e * k_E - g * k_M)/ (e + g); % 1/d, specific growth rate
      p_C = E * (k_E - r);              % J/d, mobilisation rate
      
      dE = f * p_Am * V - p_C;          % J/d, change in reserve
      dL = r * L/ 3;                    % cm/d, change in length
      dE_R = (1 - kap) * p_C - p_J;     % J/d, change in reprod buffer
      
      dER = dE_R/ V - r * E_R/ V;       % J/d.cm^3, change in [E_R]

      daELE = [1; dE; dL; dE_R]/ dER;   
    end
      
    function daVEL = dget_aVEL(H, aVEL, k_J, k_E, v, E_m, L_m, kap, g, E_G, kap_V)
      % pupal development
      % H: maturity of imago (J)
       
      % unpack variables
      a  = aVEL(1); % d, time since pupation
      V  = aVEL(2); % cm^3, structural volume of larva
      E  = aVEL(3); % J, reserve of larva
      L  = aVEL(4); % cm, structural length of imago
            
      dV = - V * k_E;                % cm^3/d, change larval structure
      e = E/ L^3/ E_m;               % -, scaled reserve density
      r = v * (e/ L - 1/ L_m)/ (e + g); % 1/d, specific growth rate
      p_C = E * (v/ L - r);          % J/d, mobilisation rate
      dE = dV * E_G * kap_V - p_C; % J/d, change in reserve
      dL = r * L/ 3;                 % cm/d, change in length
      dH = (1 - kap) * p_C - k_J * H;% J/d, change in maturity
      
      daVEL = [1; dV; dE; dL]/ dH;
    end

    function dELE = dget_ELE_W(W, ELE, f, k_M, k_E, p_J, p_Am, E_m, g, kap, L_3, f_4, w_E, mu_E, d_V)
      % larval development
      % W: wet weight of larva 
      
      % unpack variables
      E  = ELE(1); % J, reserve
      L  = ELE(2); % cm, structural length
      E_R= ELE(3); % J, reproduction buffer
      
      if L > L_3 % less feeding in instar 4
        f = f_4;
      end
      
      V = L^3;                          % cm^3, structural volume
      e = E/ V/ E_m;                    % -, scaled reserve density
      r = (e * k_E - g * k_M)/ (e + g); % 1/d, specific growth rate
      p_C = E * (k_E - r);              % J/d, mobilisation rate
      
      dE = f * p_Am * V - p_C;          % J/d, change in reserve
      dL = r * L/ 3;                    % cm/d, change in length
      dE_R = (1 - kap) * p_C - p_J;     % J/d, change in reprod buffer

      dW = 3 * L^2 * dL + (dE + dE_R) * w_E/ mu_E/ d_V; % g/d, change in wet weight
      
      dELE = [dE; dL; dE_R]/ dW;   
    end

    function dHVEL = dget_HVEL_W(W, HVEL, k_J, k_E, v, E_m, L_m, kap, g, mu_E, kap_V, w_E, d_E)
      % pupal development
      % W: wet weight of pupa (g)
       
      % unpack variables
      H  = HVEL(1); % J, maturity
      V  = HVEL(2); % cm^3, structural volume of larva
      E  = HVEL(3); % J, reserve of larva
      L  = HVEL(4); % cm, structural length of imago
      
      E_G = g * kap * E_m;           % J/cm^3, cost for structure
      dV = - V * k_E;                % cm^3/d, change larval structure
      e = E/ L^3/ E_m;               % -, scaled reserve density
      r = v * (e/ L - 1/ L_m)/ (e + g); % 1/d, specific growth rate
      p_C = E * (v/ L - r);          % J/d, mobilisation rate
      dE = dV * E_G * kap_V - p_C; % J/d, change in reserve
      dL = r * L/ 3;                 % cm/d, change in length
      dH = (1 - kap) * p_C - k_J * H;% J/d, change in maturity
      
      dW = 3 * L^2 * dL + dV + dE * w_E/ mu_E/ d_E; % g/d, change in wet weight
     
      dHVEL = [dH; dV; dE; dL]/ dW;
    end
