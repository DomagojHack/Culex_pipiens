function [data, auxData, metaData, txtData, weights] = mydata_Culex_pipiens

% File based on existing Aedes aegypti pet in AmP


%% set metaData
metaData.phylum     = 'Arthropoda'; 
metaData.class      = 'Insecta'; 
metaData.order      = 'Diptera'; 
metaData.family     = 'Culicidae';
metaData.species    = 'Culex_pipiens'; 
metaData.species_en = 'Common house mosquito'; 
metaData.ecoCode.climate = {'C'};
metaData.ecoCode.ecozone = {'TH'};
metaData.ecoCode.habitat = {'0eFm', 'Fp', 'eiTf', 'eiTi'}; %% what about urban habitats??
metaData.ecoCode.embryo  = {'Fpf'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'bjD', 'bjCi', 'eiTv'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(28); % K, body temp ??? how to estimte fro cold blooded organisms?
metaData.data_0     = {'ab'; 'aj'; 'ae'; 'am'; 'E0'; 'L_t'; 'Ww0'; 'Wwj'; 'Wwe'; 'Ri'}; 
metaData.data_1     = {'t-Ww'; 'Ww-JO'}; 


