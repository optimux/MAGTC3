%This function is used toread information from output of alphaMELTS
%The most important info are Temperature, mass, conposition, water of
%OL,CPX,OPX,PL,ILM and melt; while viscosity and Mg# of melt also counts.

%Created 2020-4-1

%NOTE: Phase_main_tbl.txt has all necessary info in magma evolution.

rootpath='F:\HT\Kuritani\';
pickfile=[rootpath,'Phase_main_tbl.txt'];

cmpdoc='E:\MELTS\CMPpick.txt';                                             %mineral composition pickup file used by column_pick.command
cmpdelim='Delimiter: space';
cmphead='Header: Matlab';

%% ========================= LIQUID INFO ==================================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: liquid\r\n');
fprintf(fid,'Columns: Temperature mass viscosity SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O Mg#\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);%read row by row
tempdata = tempdata{:}';

L=length(tempdata);

melt=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'mu',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1),...
    'NMg',zeros(L-1,1));

OL=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1));

OPX=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1));

CPX=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1));

PL=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1));

ILM=struct('T',zeros(L-1,1),...
    'mass',zeros(L-1,1),...
    'SiO2',zeros(L-1,1),...
    'TiO2',zeros(L-1,1),...
    'Al2O3',zeros(L-1,1),...
    'Fe2O3',zeros(L-1,1),...
    'FeO',zeros(L-1,1),...
    'MnO',zeros(L-1,1),...
    'MgO',zeros(L-1,1),...
    'CaO',zeros(L-1,1),...
    'Na2O',zeros(L-1,1),...
    'K2O',zeros(L-1,1),...
    'P2O5',zeros(L-1,1),...
    'H2O',zeros(L-1,1));


for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
melt.T(i-1)=str2double(mindata{1});                                        %melt Temperature [K]
melt.mass(i-1)=str2double(mindata{2});                                     %melt mass [g]
melt.mu(i-1)=str2double(mindata{3});                                       %melt viscosity [Pa.sec]

melt.SiO2(i-1)=str2double(mindata{4});                                     %SiO2 wt%
melt.TiO2(i-1)=str2double(mindata{5});                                     %TiO2
melt.Al2O3(i-1)=str2double(mindata{6});                                    %Al2O3
melt.Fe2O3(i-1)=str2double(mindata{7});                                    %Fe2O3
melt.FeO(i-1)=str2double(mindata{8});                                      %FeO
melt.MnO(i-1)=str2double(mindata{9});                                      %MnO
melt.MgO(i-1)=str2double(mindata{10});                                     %MgO
melt.CaO(i-1)=str2double(mindata{11});                                     %CaO
melt.Na2O(i-1)=str2double(mindata{12});                                    %Na2O
melt.K2O(i-1)=str2double(mindata{13});                                     %K2O
melt.P2O5(i-1)=str2double(mindata{14});                                    %P2O5
melt.H2O(i-1)=str2double(mindata{15});                                     %H2O
melt.NMg(i-1)=str2double(mindata{16});                                     %NMg
end

%% ========================= OLIVINE INFO =================================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: olivine\r\n');
fprintf(fid,'Columns: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);
tempdata = tempdata{:}';

L=length(tempdata);

for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
OL.T(i-1)=str2double(mindata{1});                                          %Temperature [K]
OL.mass(i-1)=str2double(mindata{2});                                       %mass [g]

OL.SiO2(i-1)=str2double(mindata{3});                                       %SiO2 wt%
OL.TiO2(i-1)=str2double(mindata{4});                                       %TiO2
OL.Al2O3(i-1)=str2double(mindata{5});                                      %Al2O3
OL.Fe2O3(i-1)=str2double(mindata{6});                                      %Fe2O3
OL.FeO(i-1)=str2double(mindata{7});                                        %FeO
OL.MnO(i-1)=str2double(mindata{8});                                        %MnO
OL.MgO(i-1)=str2double(mindata{9});                                        %MgO
OL.CaO(i-1)=str2double(mindata{10});                                       %CaO
OL.Na2O(i-1)=str2double(mindata{11});                                      %Na2O
OL.K2O(i-1)=str2double(mindata{12});                                       %K2O
OL.P2O5(i-1)=str2double(mindata{13});                                      %P2O5
OL.H2O(i-1)=str2double(mindata{14});                                       %H2O
end

%% ======================= ORTHOPYROXENE INFO =============================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: orthopyroxene\r\n');
fprintf(fid,'Columns: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);
tempdata = tempdata{:}';

L=length(tempdata);

for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
OPX.T(i-1)=str2double(mindata{1});                                         %Temperature [K]
OPX.mass(i-1)=str2double(mindata{2});                                      %mass [g]

OPX.SiO2(i-1)=str2double(mindata{3});                                      %SiO2 wt%
OPX.TiO2(i-1)=str2double(mindata{4});                                      %TiO2
OPX.Al2O3(i-1)=str2double(mindata{5});                                     %Al2O3
OPX.Fe2O3(i-1)=str2double(mindata{6});                                     %Fe2O3
OPX.FeO(i-1)=str2double(mindata{7});                                       %FeO
OPX.MnO(i-1)=str2double(mindata{8});                                       %MnO
OPX.MgO(i-1)=str2double(mindata{9});                                       %MgO
OPX.CaO(i-1)=str2double(mindata{10});                                      %CaO
OPX.Na2O(i-1)=str2double(mindata{11});                                     %Na2O
OPX.K2O(i-1)=str2double(mindata{12});                                      %K2O
OPX.P2O5(i-1)=str2double(mindata{13});                                     %P2O5
OPX.H2O(i-1)=str2double(mindata{14});                                      %H2O
end

%% ======================= CLINOPYROXENE INFO =============================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: clinopyroxene\r\n');
fprintf(fid,'Columns: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);
tempdata = tempdata{:}';

L=length(tempdata);

for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
CPX.T(i-1)=str2double(mindata{1});                                         %Temperature [K]
CPX.mass(i-1)=str2double(mindata{2});                                      %mass [g]

CPX.SiO2(i-1)=str2double(mindata{3});                                      %SiO2 wt%
CPX.TiO2(i-1)=str2double(mindata{4});                                      %TiO2
CPX.Al2O3(i-1)=str2double(mindata{5});                                     %Al2O3
CPX.Fe2O3(i-1)=str2double(mindata{6});                                     %Fe2O3
CPX.FeO(i-1)=str2double(mindata{7});                                       %FeO
CPX.MnO(i-1)=str2double(mindata{8});                                       %MnO
CPX.MgO(i-1)=str2double(mindata{9});                                       %MgO
CPX.CaO(i-1)=str2double(mindata{10});                                      %CaO
CPX.Na2O(i-1)=str2double(mindata{11});                                     %Na2O
CPX.K2O(i-1)=str2double(mindata{12});                                      %K2O
CPX.P2O5(i-1)=str2double(mindata{13});                                     %P2O5
CPX.H2O(i-1)=str2double(mindata{14});                                      %H2O
end

%% ========================= PLAGIOCLASE INFO =============================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: feldspar\r\n');
fprintf(fid,'Columns: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);
tempdata = tempdata{:}';

L=length(tempdata);

for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
PL.T(i-1)=str2double(mindata{1});                                          %Temperature [K]
PL.mass(i-1)=str2double(mindata{2});                                       %mass [g]

PL.SiO2(i-1)=str2double(mindata{3});                                       %SiO2 wt%
PL.TiO2(i-1)=str2double(mindata{4});                                       %TiO2
PL.Al2O3(i-1)=str2double(mindata{5});                                      %Al2O3
PL.Fe2O3(i-1)=str2double(mindata{6});                                      %Fe2O3
PL.FeO(i-1)=str2double(mindata{7});                                        %FeO
PL.MnO(i-1)=str2double(mindata{8});                                        %MnO
PL.MgO(i-1)=str2double(mindata{9});                                        %MgO
PL.CaO(i-1)=str2double(mindata{10});                                       %CaO
PL.Na2O(i-1)=str2double(mindata{11});                                      %Na2O
PL.K2O(i-1)=str2double(mindata{12});                                       %K2O
PL.P2O5(i-1)=str2double(mindata{13});                                      %P2O5
PL.H2O(i-1)=str2double(mindata{14});                                       %H2O
end

%% ========================== ILMENITE INFO ===============================
fid=fopen(cmpdoc,'w');
if(fid<3);fprintf('Failed to open CMP pick file!\n');return;end
fprintf(fid,'%s\r\n',cmpdelim);
fprintf(fid,'%s\r\n',cmphead);
fprintf(fid,'\r\n');
fprintf(fid,'File: %s\r\n',pickfile);
fprintf(fid,'Table: spinel\r\n');
fprintf(fid,'Columns: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O\r\n');
fclose(fid);

tempdata = perl('column_pick.command',cmpdoc,'> ');

if (ispc); delim = '\r\n'; else delim = '\n'; end;

tempdata = textscan(tempdata, '%s', 'Delimiter', delim);
tempdata = tempdata{:}';

L=length(tempdata);

for i=2:L
mindata=textscan(tempdata{i},'%s','Delimiter',' ');
mindata=mindata{:}';
ILM.T(i-1)=str2double(mindata{1});                                         %Temperature [K]
ILM.mass(i-1)=str2double(mindata{2});                                      %mass [g]

ILM.SiO2(i-1)=str2double(mindata{3});                                      %SiO2 wt%
ILM.TiO2(i-1)=str2double(mindata{4});                                      %TiO2
ILM.Al2O3(i-1)=str2double(mindata{5});                                     %Al2O3
ILM.Fe2O3(i-1)=str2double(mindata{6});                                     %Fe2O3
ILM.FeO(i-1)=str2double(mindata{7});                                       %FeO
ILM.MnO(i-1)=str2double(mindata{8});                                       %MnO
ILM.MgO(i-1)=str2double(mindata{9});                                       %MgO
ILM.CaO(i-1)=str2double(mindata{10});                                      %CaO
ILM.Na2O(i-1)=str2double(mindata{11});                                     %Na2O
ILM.K2O(i-1)=str2double(mindata{12});                                      %K2O
ILM.P2O5(i-1)=str2double(mindata{13});                                     %P2O5
ILM.H2O(i-1)=str2double(mindata{14});                                      %H2O
end

save([rootpath,'melt.mat'],'-struct','melt','-ascii');                     %struct name must be string, so like 'melt'
save([rootpath,'OL.mat'],'-struct','OL','-ascii');
save([rootpath,'OPX.mat'],'-struct','OPX','-ascii');
save([rootpath,'CPX.mat'],'-struct','CPX','-ascii');
save([rootpath,'PL.mat'],'-struct','PL','-ascii');
save([rootpath,'ILM.mat'],'-struct','ILM','-ascii');