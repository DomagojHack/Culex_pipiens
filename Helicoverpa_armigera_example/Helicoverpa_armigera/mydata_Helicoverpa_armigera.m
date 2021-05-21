function [data, auxData, metaData, txtData, weights] = mydata_Helicoverpa_armigera 

%% set metaData
metaData.phylum     = 'Arthropoda'; 
metaData.class      = 'Insecta'; 
metaData.order      = 'Lepidoptera'; 
metaData.family     = 'Noctuidae';
metaData.species    = 'Helicoverpa_armigera'; 
metaData.species_en = 'Cotton Bollworm'; 

metaData.ecoCode.climate = {'A','BWh','Cfa','Dfb'}; 
metaData.ecoCode.ecozone = {'TN','TH','TPa','TA'}; 
metaData.ecoCode.habitat = {'0iTg'}; 
metaData.ecoCode.embryo  = {'Tt'};
metaData.ecoCode.migrate = {'Ml'}; 
metaData.ecoCode.food    = {'bjHl','eiHn'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};  

metaData.T_typical  = C2K(25); % K, 
metaData.data_0     = {'tj'; 'te'; 'am'; 'Lb'; 'Wwj'; 'Ni'};
metaData.data_1     = {'T-ab'; 'T-tj'; 'T-te'; 't-L'}; 

metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Andre Gergs','Christian Baden'};    
metaData.date_subm = [2020 10 02];              
metaData.email    = {'andre.gergs@bayer.com'};            
metaData.address  = {'Bayer AG'};   

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2020 10 02]; 

%% set data
% zero-variate data 
data.tj = 12.7; units.tj = 'd'; label.tj = 'time since birth at pupation'; bibkey.tj = 'BarboVali2016';   
  temp.tj = C2K(27);  units.temp.tj = 'K'; label.temp.tj = 'temperature';  
data.te = 12.35; units.te = 'd'; label.te = 'duration of pupa stage'; bibkey.te = 'MiroSavo2008';   
  temp.te = C2K(25); units.temp.te = 'K'; label.temp.te = 'temperature';
data.am = 23.11; units.am = 'd'; label.am = 'life span as female imago'; bibkey.am = 'MiroSavo2008';   
  temp.am = C2K(25); units.temp.am = 'K'; label.temp.am = 'temperature';  

data.Lb = 0.16; units.Lb = 'cm'; label.Lb = 'length at birth'; bibkey.Lb = 'PatiKane2018'; 
 
data.Wwj = 0.3581;   units.Wwj = 'g'; label.Wwj = 'wet weight of pupa';      bibkey.Wwj = 'BarboVali2016';
  comment.Wwj = 'weight measure for artificial diet';  
 
data.Ni = 656.27; units.Ni = '#'; label.Ni = 'number of eggs'; bibkey.Ni = 'MiroSavo2008';   
    comment.Ni = 'average egg number excluding outliers';
 
% uni-variate data
% time - age at birth 
data.Tab = [ ... % temperature (C), age at birth (d)
12.5 19.11
15   13.57
17.5 6.6
20   5.12
25   3.31
27.5 2.4
30   2.1
32.5 2.01
35   2.04
37.5 2.11];
units.Tab  = {'C', 'd'}; label.Tab = {'temperature', 'age at birth'};  
bibkey.Tab = 'MiroSavo2008';

% temperature - duration of instars
data.TtI1 = [ ... % temperature (C), duration instar 1 (d)
15   9.31
17.5 4.3
20   4.17
25   3.27
27.5 2.46
30   2.22
32.5 2.04
35   2
37.5 2.01];
units.TtI1  = {'C', 'd'}; label.TtI1 = {'temperature', 'duration instar 1'};  
bibkey.TtI1 = 'MiroSavo2008';
%
data.TtI2 = [ ... % temperature (C), duration instar 2 (d)
15   7.87 
17.5 3.14
20   3.06
25   1.91
27.5 1.91
30   1.6
32.5 1.21
35   1.02
37.5 1.19];
units.TtI2  = {'C', 'd'}; label.TtI2 = {'temperature', 'duration instar 2'};  
bibkey.TtI2 = 'MiroSavo2008';
%
data.TtI3 = [ ... % temperature (C), duration instar 3 (d)
15   8.22
17.5 3.84
20   2.98
25   2.03
27.5 2.03
30   1.96
32.5 1.57
35   1.12
37.5 1.76];
units.TtI3  = {'C', 'd'}; label.TtI3 = {'temperature', 'duration instar 3'};  
bibkey.TtI3 = 'MiroSavo2008';
%
data.TtI4 = [ ... % temperature (C), duration instar 4 (d)
15   10.41
17.5 5.14
20   4.06
25   2.31
27.5 3.53
30   3.03
32.5 2.55
35   2.49
37.5 3.05];
units.TtI4  = {'C', 'd'}; label.TtI4 = {'temperature', 'duration instar 4'};  
bibkey.TtI4 = 'MiroSavo2008';
%
data.TtI5 = [ ... % temperature (C), duration instar 5 (d)
15   19.81
17.5 10.74
20   8.6
25   5.18
27.5 5.06
30   4
32.5 3.37
35   4.51
37.5 3.75];
units.TtI5  = {'C', 'd'}; label.TtI5 = {'temperature', 'duration instar 5'};  
bibkey.TtI5 = 'MiroSavo2008';

% temperature - duration of pupa stage
data.Tte = [ ... % temperature (C), duration pupa stage (d)
17.5 23.81
20   20.8
25   12.35
27.5 11.29
30   9.07
32.5 8.36
35   8.56
37.5 8.33];
units.Tte  = {'C', 'd'}; label.Tte = {'temperature', 'duration pupa stage'};  
bibkey.Tte = 'MiroSavo2008';

% time - length
data.tL = [ ... % time (d), length (cm)
1 0.0263
2 0.0304
3 0.0386
4 0.0421
5 0.062
7 0.137
8 0.154
9 0.197
10 0.249];
units.tL   = {'d', 'cm'};  label.tL = {'time since birth', 'head length'};  
temp.tL    = C2K(26);  units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = 'BarboVali2016';

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.v = 3 * weights.psd.v;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Discussion points

%% Links
metaData.links.id_CoL = '261649250049fd8310aec4189043fd3b'; % Cat of Life
metaData.links.id_EoL = '2869488'; % Ency of Life
metaData.links.id_Wiki = 'Helicoverpa_armigera'; % Wikipedia
metaData.links.id_ADW = ''; % ADW, not present
metaData.links.id_Taxo = '81191'; % Taxonomicon

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Helicoverpa_armigera}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'BarboVali2016'; type = 'Article'; bib = [ ... 
'author = {Barboisa, T. A. M. and Mendes, S. M. and Rodrigues, G. T. and Ribeiro, P. E. d. A. and Santos, C. A. d. and Valicente, F. H. and Oliveira, C. M. d.}, ' ...
'year = {2016}, ' ...
'title = {Comparison of biology between \emph{Helicoverpa zea} and \emph{Helicoverpa armigera} ({L}epidoptera: {N}octuidae) reared on artificial diets}, ' ...
'journal = {Florida Entomologist}, ' ...
'volume = {99}, ' ...
'pages = {72-76}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PatiKane2018'; type = 'Article'; bib = [ ... 
'author = {Patil, V. M. and Patel, Z. P. and Oak, P. S. and Chauhan, R. C. and Kaneriya, P. B.}, ' ...
'year = {2018}, ' ...
'title = {Biology of fruit borer, \emph{Helicoverpa armigera} ({H}ubner) in/on chilli fruits}, ' ...
'journal = {International Journal of Entomology Research}, ' ...
'volume = {3}, ' ...
'pages = {6-12}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'MiroSavo2008'; type = 'Article'; bib = [ ... 
'author = {Mironidis, G. K. and Savopoulou-Soultani, M.}, ' ...
'year = {2008}, ' ...
'title = {Development, Survivorship, and Reproduction of \emph{Helicoverpa armigera} ({L}epidoptera: {N}octuidae) Under Constant and Alternating Temperatures}, ' ...
'journal = {Physiological Ecology}, ' ...
'volume = {104}, ' ...
'pages = {751-764}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

