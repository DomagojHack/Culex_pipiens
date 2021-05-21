function [data, auxData, metaData, txtData, weights] = mydata_Culex_pipiens

% File based on existing Aedes aegypti and Helicoverpa armigera pet in AmP 


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


% 0	start of development	j	end of acceleration	R	first breeding
% h	hatching (leaving egg)	p	puberty (start of allocation to reproduction)	e	emergence (imago)
% b	birth (start of feeding)	i	ultimate state	s	start of acceleration
% x	weaning (or fledging)	m	death (life span

% ab age at birth
% aj age at end of acceleration
% ae age at emergance
% am age at death
% E0 reserve energy at 0
% L_t lenght at time t
% Ww0 weight at 0
% WWj weight at j
% Wwe weight at e
% Ri reproduction rate at ultimate state (i)

% t-Ww timem wet weight
% Ww-JO Wet weight O2 consumtion

data.ab = 27.26; units.ab = 'h'; label.ab = 'age at birth';   bibkey.ab = "Subra1980"
temp.ab = C2K(28.1);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.tj = 8;     units.tj = 'd';     label.tj = 'time since birth at pupation'; bibkey.tj = 'Subra1980';   
temp.tj = C2K(28);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
data.te = 1.666;     units.te = 'd';     label.te = 'time since pupation at emergence'; bibkey.te = 'Subra1980';
temp.te = C2K(28);  units.temp.te = 'K'; label.temp.te = 'temperature';
data.am = 22.8;    units.am = 'd';     label.am = 'life span as imago';       bibkey.am = 'Andreadis2014';   
temp.am = C2K(27.5);  units.temp.am = 'K'; label.temp.am = 'temperature';
data.t1 = 1.2;    units.t1 = 'd';     label.t1 = 'duration of instar 1';         bibkey.t1 = 'Loetti2011';   
  temp.t1 = C2K(30);  units.temp.t1 = 'K'; label.temp.t1 = 'temperature'; 
data.t2 = 1.1;    units.t2 = 'd';     label.t2 = 'duration of instar 2';         bibkey.t2 = 'Loetti2011';   
  temp.t2 = C2K(30);  units.temp.t2 = 'K'; label.temp.t2 = 'temperature'; 
data.t3 = 1.2;    units.t3 = 'd';     label.t3 = 'duration of instar 3';         bibkey.t3 = 'Loetti2011';   
  temp.t3 = C2K(30);  units.temp.t3 = 'K'; label.temp.t3 = 'temperature'; 
data.t4 = 2.7;    units.t4 = 'd';     label.t4 = 'duration of instar 4';         bibkey.t4 = 'Loetti2011';   
  temp.t4 = C2K(30);  units.temp.t4 = 'K'; label.temp.t4 = 'temperature'; 

% Still from Aedes egypti
data.L1 = 0.278;  units.L1 = 'mm'; label.L1 = 'length head capsule of instar 1'; bibkey.L1 = 'Chri1960';
data.L2 = 0.463;  units.L2 = 'mm'; label.L2 = 'length head capsule of instar 2'; bibkey.L2 = 'Chri1960';
data.L3 = 0.74;   units.L3 = 'mm'; label.L3 = 'length head capsule of instar 3'; bibkey.L3 = 'Chri1960';
data.L4 = 0.98;   units.L4 = 'mm'; label.L4 = 'length head capsule of instar 4'; bibkey.L4 = 'Chri1960';

% Still from Aedes egypti
data.Ww0 = 0.0137; units.Ww0 = 'mg'; label.Ww0 = 'initial wet weight';      bibkey.Ww0 = 'Chri1960';
  comment.Ww0 = 'actually 0.01 in Chris1960 but this is less than L1';
data.Wwj = 4.752;  units.Wwj = 'mg'; label.Wwj = 'wet weight of pupa';      bibkey.Wwj = 'Chri1960';
data.Wwe = 3.04;   units.Wwe = 'mg'; label.Wwe = 'wet weight of imago';     bibkey.Wwe = 'Chri1960';
data.E0 = 18.2*4.184; units.E0 = 'mJ'; label.E0 = 'initial energy content'; bibkey.E0 = 'Brie1990';
data.Ri  = 400/21; units.Ri  = '#/d';  label.Ri  = 'maximum reprod rate';   bibkey.Ri  = 'Chri1960';   
  temp.Ri = C2K(28); units.temp.Ri = 'K'; label.temp.Ri = 'temperature';
 
% uni-variate data

% duration of L1 instalr and temperature

data.TtI1= [
  7   6.3
  10  4.7
  15  4.2
  20  2.8
  25  2.0
  30  1.2
  33  1.3
];
units.TtI1 = {'C', 'd'}; label.TtI1 = {"temperature", "duration instar 1"};
bibkey.TtI1 = "Loetti2011";

% duration of L2 instalr and temperature

data.TtI2= [
  7   12.7
  10  8.6
  15  3.2
  20  1.7
  25  1.5
  30  1.1
  33  1.3
];
units.TtI2 = {'C', 'd'}; label.TtI2 = {"temperature", "duration instar 2"};
bibkey.TtI2 = "Loetti2011";

% duration of L3 instalr and temperature

data.TtI3= [
  7   10.1
  10  7.1
  15  3.2
  20  2.1
  25  1.7
  30  1.2
  33  2.3
];
units.TtI3 = {'C', 'd'}; label.TtI3 = {"temperature", "duration instar 3"};
bibkey.TtI3 = "Loetti2011";


% duration of L4 instalr and temperature

data.TtI4= [
  7   10.1
  10  7.1
  15  3.2
  20  2.1
  25  1.7
  30  1.2
  33  2.3
];
units.TtI4 = {'C', 'd'}; label.TtI4 = {"temperature", "duration instar 4"};
bibkey.TtI4 = "Loetti2011";


% duration of pupa and temperature

% data.TtP= [
%   7   10.1
%   10  7.1
%   15  3.2
%   20  2.1
%   25  1.7
%   30  1.2
%   33  2.3
% ];
% units.TtP = {'C', 'd'}; label.TtP = {"temperature", "duration pupa"};
% bibkey.TtP = "Loetti2011";


% duration of larva to adult emergence and temperature

% data.Tte= [
%   10  39.8
%   15  21.6
%   20  13.3
%   25  10.2
%   30  8.0
%   33  11.5
% ];
% units.Tte = {'C', 'd'}; label.Tte = {"temperature", "duration larvae 1 to adult"};
% bibkey.Tte = "Loetti2011";




% % time - wet weights of larva  
% data.tW =  [ ... % time since birth (d), wet weight (mg)
%         0           0.0137
%         0.75       	0.065
%         0.958333333	0.089
%         1           0.09
%         1.041666667	0.118
%         1.083333333	0.119
%         1.166666667	0.146
%         1.833333333	0.414
%         1.875       0.54
%         1.916666667	0.52
%         2           0.55
%         2.125   	0.58
%         2.208333333	0.69
%         2.25    	0.73
%         2.875       1.71
%         2.916666667	1.76
%         2.958333333	1.5
%         3.166666667	2.07
%         3.958333333	2.96
%         4.791666667	3.59
%         4.875       3.86
%         4.958333333	4.13
%         6           5.6];
% units.tW  = {'d', 'mg'}; label.tW = {'time since birth', 'wet weight'};  
% temp.tW   = C2K(28);  units.temp.tW = 'K'; label.temp.tW = 'temperature';
% bibkey.tW = 'Chri1960';

% % time - wing-lenght 

% % respiration data		
% TWJO = [ ... % temp (C), wet weight (mg), O2 consumption mm^3/g/h 
%  25.5  10.1   990   % larva
%  27.9  40 	 1200   % larva
%  27.6  45.8  1048   % larva
%  27.5  14.8   736   % pupa
%  28    31.7   902   % pupa
%  28     4.86 4713   % imago
%  28 	6.38 3495]; % imago
% data.WJO = TWJO(:,[2 3]);
% units.WJO   = {'mg', 'mm^3/h.g'};  label.WJO = {'wet weight', 'O_2 consumption'};  
% temp.WJO    = C2K(TWJO(:,1));  units.temp.WJO = 'K'; label.temp.WJO = 'temperature';
% bibkey.WJO = 'Chri1960';
  
% %% set weights for all real data
weights = setweights(data, []);
% weights.tW = 5 * weights.tW;
% weights.Ww0 = 10 * weights.Ww0;
% weights.Wwj = 5 * weights.Wwj;

% %% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.v = 5 * weights.psd.v;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

% %% Discussion points
% D1 = 'Instar 4 grows slower than 1-3, probably due to ceasing of feeding in preparation for pupation, modelled as lower f';
% metaData.discussion = struct('D1', D1);

% %% Facts
% F1 = '4 instars between birth and pupation';
% metaData.bibkey.F1 = 'Chri1960'; 
% F2 = 'female emago takes blood for eggs, but larvae also allocate to reproduction';
% metaData.bibkey.F2 = 'Chri1960'; 
% metaData.facts = struct('F1',F1,'F2',F2);

% % %% Links
% % metaData.links.id_CoL = '659d82334ce06794ac14a699fe41bb4d'; % Cat of Life
% % metaData.links.id_EoL = '41592971'; % Ency of Life
% % metaData.links.id_Wiki = 'Aedes_aegypti'; % Wikipedia
% % metaData.links.id_ADW = 'Aedes_aegypti'; % ADW
% % metaData.links.id_Taxo = '28492'; % Taxonomicon

% %% References
% bibkey = 'Wiki'; type = 'Misc'; bib = ...
% 'howpublished = {\url{http://en.wikipedia.org/wiki/Aedes_aegypti}}';
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
% bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
% 'author = {Kooijman, S.A.L.M.}, ' ...
% 'year = {2010}, ' ...
% 'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
% 'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
% 'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
% 'howpublished = {\url{../../../bib/Kooy2010.html}}'];
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
% bibkey = 'Chri1960'; type = 'Book'; bib = [ ... 
% 'author = {Christophers, S. R.}, ' ... 
% 'year = {1960}, ' ...
% 'title = {\emph{Aedes aegypti} ({L}.) the Yellow Fever Mosquito: its Life History, Bionomics and Structure}, ' ...
% 'publisher = {Cambridge University Press}, ' ...
% 'address = {Cambridge}'];
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
% bibkey = 'Brie1990'; type = 'Article'; bib = [ ... 
% 'author = {Briegel, H.}, ' ... 
% 'year = {1990}, ' ...
% 'title = {Metabolic relationship between female body size, reserves, and fecundity of \emph{Aedes aegypti}}, ' ...
% 'journal = {J. Insect Physiol.}, ' ...
% 'volume = {36}, ' ...
% 'pages = {165--172}'];
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %
% bibkey = 'Subra1980'; type = 'Book'; bib = [ ... 
% 'author = {Christophers, S. R.}, ' ... 
% 'year = {1960}, ' ...
% 'title = {\emph{Aedes aegypti} ({L}.) the Yellow Fever Mosquito: its Life History, Bionomics and Structure}, ' ...
% 'publisher = {Cambridge University Press}, ' ...
% 'address = {Cambridge}'];
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
