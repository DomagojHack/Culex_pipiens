function [prdData, info] = predict_Helicoverpa_armigera(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  if s_1 < 1 || s_2 < 1 || s_3 < 1 || s_4 < 1 
    prdData = []; info = 0; return
  end
  
  % compute temperature correction factors
  TC_27 = tempcorr(temp.tj , T_ref, [T_A; T_L; T_AL]);           % data at 27 C
  TC_25 = tempcorr(temp.te, T_ref, [T_A; T_L; T_AL]);           % data at 25 C
  TC_ab = tempcorr(C2K(Tab(:,1)), T_ref, [T_A; T_L; T_AL]);     % instar 0, embryo
  TC_TL = tempcorr(C2K(TtI1(:,1)), T_ref, [T_A; T_L; T_AL]);  % Temperature correction for larvae
  TC_te = tempcorr(C2K(Tte(:,1)), T_ref, [T_A; T_L; T_AL]);    % Temperature correction for pupae
  
  TC = tempcorr(temp.tL, T_ref, [T_A; T_L; T_AL]); % time length data data at 26 C 
  vT = v * TC; kT_J = k_J * TC; kT_M = k_M * TC; pT_Am = p_Am * TC;
  
  % life cycle
  pars_tj = [g k v_Hb v_He s_j kap kap_V];
  [t_j, t_e, t_b, l_j, l_e, l_b, rho_j, v_Rj, u_Ee, info] = get_tj_hex(pars_tj, f);
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  E_0 = p_Am * initial_scaled_reserve(f, pars_UE0); % J, initial energy content
    
  % birth
  a_b = t_b/ k_M;                     % d, age at birth at f
  L_b = L_m * l_b;                    % cm, structural length at birth
  Lw_b = L_b/ del_M;                  % cm, head length at birth
  Ww_b =  L_b^3 * (1 + f * w);        % g, wet weight at birth

  % metam
  E_R_ref = (1 - kap) * E_m * (g + l_b)/ (1 - l_b); % J/cm^3, reference value for [E_R]
  E_Rj = s_j * E_R_ref;               % J/cm^3, threshold reprod buffer density at pupation
  tT_j = (t_j - t_b)/ k_M/ TC_27;     % d, time since birth at pupation
  L_j = L_m * l_j;                    % cm, struc length at metam
  Lw_j = L_j/ del_M;                  % cm, length at metam
  s_M = L_j/ L_b;                     % -, acceleration factor
  E = f * E_m * L_j^3;                % J, reserve at pupation
  E_j = E + E_Rj;                     % J, reserve at pupation, including reproduction buffer
  Ww_j = (L_j^3 + E_j * w_E/ mu_E/ d_E); % g, wet weight at pupation
  
  % emergence
  v_j = v * s_M;                      % cm/d, energy conductance of imago
  L_mj = L_m * s_M;                   % cm, max length of imago
  [H aVEL] = ode45(@dget_aVEL, [0; E_He], [0; L_j^3; f * E_m * L_j^3; 1e-4], [], k_J, v/ L_b, v_j, E_m, L_mj, kap, g, E_G, kap_V);
  V_e = aVEL(end,2);                  % cm^3, larval structure at emergence (should be close to zero)
  E_e = aVEL(end,3);                  % J, reserve at emergence
  L_e = aVEL(end,4);                  % cm, length at emergence
  E_e = E_e + kap_V * E_G * V_e;      % J, convert last bit of larval structure to reserve instantaneously
  tT_e = (t_e - t_j)/ k_M/ TC_25;     % d, time since pupation at emergence
  
  % instar age, length
  Lw_1 = Lw_b; Lw_2 = Lw_1 * s_1^0.5; Lw_3 = Lw_2 * s_2^0.5; Lw_4 = Lw_3 * s_3^0.5; Lw_5 = Lw_4 * s_4^0.5;  % cm, head length of instar 1-4
  rT_j = TC_25 * g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_1 = log(Lw_2/ Lw_1) * 3/ rT_j;  tT_2 = log(Lw_3/ Lw_2) * 3/ rT_j;
  tT_3 = log(Lw_4/ Lw_3) * 3/ rT_j;  tT_4 = log(Lw_5/ Lw_4) * 3/ rT_j;
  tT_5 = log(Lw_j/ Lw_5) * 3/ rT_j;  
  % reproduction 
  N_i = kap_R * L_j^3 * E_Rj/ E_0;     % #, life time reprod output

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);     % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_25;              % d, mean life span at T
  
  % pack to output
  prdData.tj = tT_j;
  prdData.te = tT_e;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Wwj = Ww_j;
  prdData.Ni = N_i;
  
  % uni-variate data 
  
  % temperature - age at birth
  EaT_b = a_b ./ TC_ab;

  % temperature - development time instar 1
  rT_j = g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_1 = log(Lw_2/ Lw_1) * 3/ rT_j; 
  ET_tI_1 = tT_1 ./ TC_TL;    
  
  % temperature - development time instar 2
  rT_j = g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_2 = log(Lw_3/ Lw_2) * 3/ rT_j;
  ET_tI_2 = tT_2 ./ TC_TL;     
  
  % temperature - development time instar 3
  rT_j = g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_3 = log(Lw_4/ Lw_3) * 3/ rT_j; 
  ET_tI_3 = tT_3 ./ TC_TL;    
  
  % temperature - development time instar 4
  rT_j = g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_4 = log(Lw_5/ Lw_4) * 3/ rT_j; 
  ET_tI_4 = tT_4 ./ TC_TL;   
  
  % temperature - development time instar 5
  rT_j = g * k_M * (f/ l_b - 1)/ (f + g);  % 1/d, specific growth rate
  tT_5 = log(Lw_j/ Lw_5) * 3/ rT_j;  
  ET_tI_5 = tT_5 ./ TC_TL;     
  
  % temperature - development time pupa
  t_e = (t_e - t_j)/ k_M;      % d, time since pupation at emergence
  ET_e = t_e ./ TC_te;   
  
  % time - length
   ELE_0 = [f * E_m * L_b^3; L_b; 0]; % initial state vector
   pT_J = kT_J * E_Hb;                % J/d, maturity maintenance from birth till pupation
   [t ELE] = ode45(@dget_ELE, tL(:,1), ELE_0, [], f, kT_M, vT/ L_b, pT_J, pT_Am/ L_b, E_m, g, kap);
   ELw = ELE(:,2)./ del_Mh; % cm, length 

  
  %pack to output  
  prdData.Tab   = EaT_b;
  prdData.TtI1  = ET_tI_1; 
  prdData.TtI2  = ET_tI_2;
  prdData.TtI3  = ET_tI_3;
  prdData.TtI4  = ET_tI_4;
  prdData.TtI5  = ET_tI_5;
  prdData.Tte   = ET_e;
  prdData.tL  = ELw;

end

%% subfunctions    
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
    
function dELE = dget_ELE(t, ELE, f, k_M, k_E, p_J, p_Am, E_m, g, kap)
      % larval development
      % t: time since birth
      % p_Am: [p_Am] = {p_Am}/ L_b
      
      % unpack variables
      E  = ELE(1); % J, reserve
      L  = ELE(2); % cm, structural length
      %E_R = ELE(3); % J, reproduction buffer
      
      V = L^3;                          % cm^3, structural volume
      e = E/ V/ E_m;                    % -, scaled reserve density
      r = (e * k_E - g * k_M)/ (e + g); % 1/d, specific growth rate
      p_C = E * (k_E - r);             % J/d, mobilisation rate
      
      dE = f * p_Am * V - p_C;          % J/d, change in reserve
      dL = r * L/ 3;                    % cm/d, change in length
      dE_R = (1 - kap) * p_C - p_J;     % J/d, change in reprod buffer

      dELE = [dE; dL; dE_R];   
end
    

