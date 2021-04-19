%This script calculates all mass fraction of oxides in all mineral endmembers
%Created on 2020-7-7

global MOxide
global MMins
global MCSEM   %MCS in End-Members

%End-members of minerals, from Schutt&Lesher2006 and densityfull.m
%% ==========================End-member abbreviation=======================
%Abbreviations of End Members

%Mineral                   Abbr.                    Formaula
%--------------------------olivine----------------------------------
%Forsterite                Fo                       Mg2SiO4
%Fayalite                  Fa                       Fe2SiO4

%--------------------------OPX--------------------------------------
%Acmite                    OAm                     NaFe(III)Si2O6  
%Jadeite                   OJa                     NaAlSi2O6  
%Fassaite                  OFs                     CaAl2SiO6
%(Mg,Fe)-Fassaite          OMFFs                   (Mg,Fe)Fe(III)2SiO6
%Ca-Tschermak              OCaTm                   CaAl2SiO6 
%CaTi-Tschermak            OCTTm                   CaTiAl2O6
%CrCa-Tschermak            OCCTm                   CaCr2SiO6
%Mg-Tschermak              OMTm                    MgAl2SiO6
%(Mg,Fe)Ti-Tschermak       OMFTTm                  (Mg,Fe)TiAl2O6
%Fe-Tschermak              OFTm                    Fe(II)Al2SiO6
%Cr(Mg,Fe)-Tschermak       OCMFTm                  (Mg,Fe)Cr2SiO6
%Diopside                  ODi                     CaMgSi2O6  
%Hedenbergite              OHd                     CaFeSi2O6 
%Enstatite                 OEn                     Mg2Si2O6 
%Ferrosilite               OFes                    Fe2Si2O6  

%--------------------------CPX-------------------------------------
%Acmite                    CAm                     NaFe(III)Si2O6  
%Jadeite                   CJa                     NaAlSi2O6  
%Fassaite                  CFs                     CaAl2SiO6
%(Mg,Fe)-Fassaite          CMFFs                   (Mg,Fe)Fe(III)2SiO6
%Ca-Tschermak              CCaTm                   CaAl2SiO6 
%CaTi-Tschermak            CCTTm                   CaTiAl2O6
%CrCa-Tschermak            CCCTm                   CaCr2SiO6
%Mg-Tschermak              CMTm                    MgAl2SiO6
%(Mg,Fe)Ti-Tschermak       CMFTTm                  (Mg,Fe)TiAl2O6
%Fe-Tschermak              CFTm                    Fe(II)Al2SiO6
%Cr(Mg,Fe)-Tschermak       CCMFTm                  (Mg,Fe)Cr2SiO6
%Diopside                  CDi                     CaMgSi2O6  
%Hedenbergite              CHd                     CaFeSi2O6 
%Enstatite                 CEn                     Mg2Si2O6 
%Ferrosilite               CFes                    Fe2Si2O6  

%-------------------------garnet------------------------------------
%Pyrope                    Py                      Mg3Al2[SiO4]3
%Almandine                 Amd                     Fe3Al2[SiO4]3
%Spessartite               Sa                      Mn3Al2[SiO4]3
%Andradite                 Ad                      Ca3Fe(III)2[SiO4]3
%Grossular                 Gs                      Ca3Al2[SiO4]3
%Uvarovite                 Uv                      Ca3Cr2[SiO4]3

%--------------------------spinel------------------------------------
%Spinel                    Sp                      MgAl2O4
%Chromite                  Cm                      Fe(II)Cr2O4
%Hercynite                 Hy                      Fe(II)Al2O4
%Mg-chromite               MCm                     MgCr2O4
%Magnetite                 Mag                     Fe(III)Fe(II)Fe(III)O4
%Qandilite                 Qd                      Mg2TiO4  
%Ulvospinel                Us                      TiFe(II)2O4

%--------------------------feldspar----------------------------------
%Albite                    Ab                      NaAlSi3O8
%Anorthite                 An                      CaAl2Si2O8

%NOTE: in general, Fo, Fa, Di, Hd, En, Fs, Py, Gs, Sp, Us, An and Ab are most often used.

MCSEM=struct('SiO2',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'TiO2',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Al2O3',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'FeO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Fe2O3',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'MnO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'MgO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'CaO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Na2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'K2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'P2O5',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'H2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0));

%SiO2 mass fraction in all end-members
MCSEM.SiO2.Fo=MOxide.SiO2/MMins.Fo;
MCSEM.SiO2.Fa=MOxide.SiO2/MMins.Fa;
MCSEM.SiO2.Di=2.0*MOxide.SiO2/MMins.Di;
MCSEM.SiO2.Hd=2.0*MOxide.SiO2/MMins.Hd;
MCSEM.SiO2.En=2.0*MOxide.SiO2/MMins.En;
MCSEM.SiO2.Fs=2.0*MOxide.SiO2/MMins.Fs;
MCSEM.SiO2.Py=3.0*MOxide.SiO2/MMins.Py;
MCSEM.SiO2.Gs=3.0*MOxide.SiO2/MMins.Gs;
MCSEM.SiO2.Sp=0.0*MOxide.SiO2/MMins.Sp;
MCSEM.SiO2.Us=0.0*MOxide.SiO2/MMins.Us;
MCSEM.SiO2.An=2.0*MOxide.SiO2/MMins.An;
MCSEM.SiO2.Ab=3.0*MOxide.SiO2/MMins.Ab;

%TiO2 mass fraction in all end-members
MCSEM.TiO2.Fo=0.0*MOxide.TiO2/MMins.Fo;
MCSEM.TiO2.Fa=0.0*MOxide.TiO2/MMins.Fa;
MCSEM.TiO2.Di=0.0*MOxide.TiO2/MMins.Di;
MCSEM.TiO2.Hd=0.0*MOxide.TiO2/MMins.Hd;
MCSEM.TiO2.En=0.0*MOxide.TiO2/MMins.En;
MCSEM.TiO2.Fs=0.0*MOxide.TiO2/MMins.Fs;
MCSEM.TiO2.Py=0.0*MOxide.TiO2/MMins.Py;
MCSEM.TiO2.Gs=0.0*MOxide.TiO2/MMins.Gs;
MCSEM.TiO2.Sp=0.0*MOxide.TiO2/MMins.Sp;
MCSEM.TiO2.Us=MOxide.TiO2/MMins.Us;
MCSEM.TiO2.An=0.0*MOxide.TiO2/MMins.An;
MCSEM.TiO2.Ab=0.0*MOxide.TiO2/MMins.Ab;

%Al2O3 mass fraction in all end-members
MCSEM.Al2O3.Fo=0.0*MOxide.Al2O3/MMins.Fo;
MCSEM.Al2O3.Fa=0.0*MOxide.Al2O3/MMins.Fa;
MCSEM.Al2O3.Di=0.0*MOxide.Al2O3/MMins.Di;
MCSEM.Al2O3.Hd=0.0*MOxide.Al2O3/MMins.Hd;
MCSEM.Al2O3.En=0.0*MOxide.Al2O3/MMins.En;
MCSEM.Al2O3.Fs=0.0*MOxide.Al2O3/MMins.Fs;
MCSEM.Al2O3.Py=MOxide.Al2O3/MMins.Py;
MCSEM.Al2O3.Gs=MOxide.Al2O3/MMins.Gs;
MCSEM.Al2O3.Sp=MOxide.Al2O3/MMins.Sp;
MCSEM.Al2O3.Us=0.0*MOxide.Al2O3/MMins.Us;
MCSEM.Al2O3.An=MOxide.Al2O3/MMins.An;
MCSEM.Al2O3.Ab=0.5*MOxide.Al2O3/MMins.Ab;

%FeO mass fraction in all end-members
MCSEM.FeO.Fo=0.0*MOxide.FeO/MMins.Fo;
MCSEM.FeO.Fa=2.0*MOxide.FeO/MMins.Fa;
MCSEM.FeO.Di=0.0*MOxide.FeO/MMins.Di;
MCSEM.FeO.Hd=MOxide.FeO/MMins.Hd;
MCSEM.FeO.En=0.0*MOxide.FeO/MMins.En;
MCSEM.FeO.Fs=2.0*MOxide.FeO/MMins.Fs;
MCSEM.FeO.Py=0.0*MOxide.FeO/MMins.Py;
MCSEM.FeO.Gs=0.0*MOxide.FeO/MMins.Gs;
MCSEM.FeO.Sp=0.0*MOxide.FeO/MMins.Sp;
MCSEM.FeO.Us=2.0*MOxide.FeO/MMins.Us;
MCSEM.FeO.An=0.0*MOxide.FeO/MMins.An;
MCSEM.FeO.Ab=0.0*MOxide.FeO/MMins.Ab;

%Fe2O3 mass fraction in all end-members
MCSEM.Fe2O3.Fo=0.0*MOxide.Fe2O3/MMins.Fo;
MCSEM.Fe2O3.Fa=0.0*MOxide.Fe2O3/MMins.Fa;
MCSEM.Fe2O3.Di=0.0*MOxide.Fe2O3/MMins.Di;
MCSEM.Fe2O3.Hd=0.0*MOxide.Fe2O3/MMins.Hd;
MCSEM.Fe2O3.En=0.0*MOxide.Fe2O3/MMins.En;
MCSEM.Fe2O3.Fs=0.0*MOxide.Fe2O3/MMins.Fs;
MCSEM.Fe2O3.Py=0.0*MOxide.Fe2O3/MMins.Py;
MCSEM.Fe2O3.Gs=0.0*MOxide.Fe2O3/MMins.Gs;
MCSEM.Fe2O3.Sp=0.0*MOxide.Fe2O3/MMins.Sp;
MCSEM.Fe2O3.Us=0.0*MOxide.Fe2O3/MMins.Us;
MCSEM.Fe2O3.An=0.0*MOxide.Fe2O3/MMins.An;
MCSEM.Fe2O3.Ab=0.0*MOxide.Fe2O3/MMins.Ab;

%MnO mass fraction in all end-members
MCSEM.MnO.Fo=0.0*MOxide.MnO/MMins.Fo;
MCSEM.MnO.Fa=0.0*MOxide.MnO/MMins.Fa;
MCSEM.MnO.Di=0.0*MOxide.MnO/MMins.Di;
MCSEM.MnO.Hd=0.0*MOxide.MnO/MMins.Hd;
MCSEM.MnO.En=0.0*MOxide.MnO/MMins.En;
MCSEM.MnO.Fs=0.0*MOxide.MnO/MMins.Fs;
MCSEM.MnO.Py=0.0*MOxide.MnO/MMins.Py;
MCSEM.MnO.Gs=0.0*MOxide.MnO/MMins.Gs;
MCSEM.MnO.Sp=0.0*MOxide.MnO/MMins.Sp;
MCSEM.MnO.Us=0.0*MOxide.MnO/MMins.Us;
MCSEM.MnO.An=0.0*MOxide.MnO/MMins.An;
MCSEM.MnO.Ab=0.0*MOxide.MnO/MMins.Ab;

%MgO mass fraction in all end-members
MCSEM.MgO.Fo=2.0*MOxide.MgO/MMins.Fo;
MCSEM.MgO.Fa=0.0*MOxide.MgO/MMins.Fa;
MCSEM.MgO.Di=MOxide.MgO/MMins.Di;
MCSEM.MgO.Hd=0.*MOxide.MgO/MMins.Hd;
MCSEM.MgO.En=2.0*MOxide.MgO/MMins.En;
MCSEM.MgO.Fs=0.0*MOxide.MgO/MMins.Fs;
MCSEM.MgO.Py=3.0*MOxide.MgO/MMins.Py;
MCSEM.MgO.Gs=0.0*MOxide.MgO/MMins.Gs;
MCSEM.MgO.Sp=MOxide.MgO/MMins.Sp;
MCSEM.MgO.Us=0.0*MOxide.MgO/MMins.Us;
MCSEM.MgO.An=0.0*MOxide.MgO/MMins.An;
MCSEM.MgO.Ab=0.0*MOxide.MgO/MMins.Ab;

%CaO mass fraction in all end-members
MCSEM.CaO.Fo=0.0*MOxide.CaO/MMins.Fo;
MCSEM.CaO.Fa=0.0*MOxide.CaO/MMins.Fa;
MCSEM.CaO.Di=MOxide.CaO/MMins.Di;
MCSEM.CaO.Hd=MOxide.CaO/MMins.Hd;
MCSEM.CaO.En=0.0*MOxide.CaO/MMins.En;
MCSEM.CaO.Fs=0.0*MOxide.CaO/MMins.Fs;
MCSEM.CaO.Py=0.0*MOxide.CaO/MMins.Py;
MCSEM.CaO.Gs=3.0*MOxide.CaO/MMins.Gs;
MCSEM.CaO.Sp=0.0*MOxide.CaO/MMins.Sp;
MCSEM.CaO.Us=0.0*MOxide.CaO/MMins.Us;
MCSEM.CaO.An=MOxide.CaO/MMins.An;
MCSEM.CaO.Ab=0.0*MOxide.CaO/MMins.Ab;

%Na2O mass fraction in all end-members
MCSEM.Na2O.Fo=0.0*MOxide.Na2O/MMins.Fo;
MCSEM.Na2O.Fa=0.0*MOxide.Na2O/MMins.Fa;
MCSEM.Na2O.Di=0.0*MOxide.Na2O/MMins.Di;
MCSEM.Na2O.Hd=0.0*MOxide.Na2O/MMins.Hd;
MCSEM.Na2O.En=0.0*MOxide.Na2O/MMins.En;
MCSEM.Na2O.Fs=0.0*MOxide.Na2O/MMins.Fs;
MCSEM.Na2O.Py=0.0*MOxide.Na2O/MMins.Py;
MCSEM.Na2O.Gs=0.0*MOxide.Na2O/MMins.Gs;
MCSEM.Na2O.Sp=0.0*MOxide.Na2O/MMins.Sp;
MCSEM.Na2O.Us=0.0*MOxide.Na2O/MMins.Us;
MCSEM.Na2O.An=0.0*MOxide.Na2O/MMins.An;
MCSEM.Na2O.Ab=MOxide.Na2O/MMins.Ab;

%K2O mass fraction in all end-members
MCSEM.K2O.Fo=0.0*MOxide.K2O/MMins.Fo;
MCSEM.K2O.Fa=0.0*MOxide.K2O/MMins.Fa;
MCSEM.K2O.Di=0.0*MOxide.K2O/MMins.Di;
MCSEM.K2O.Hd=0.0*MOxide.K2O/MMins.Hd;
MCSEM.K2O.En=0.0*MOxide.K2O/MMins.En;
MCSEM.K2O.Fs=0.0*MOxide.K2O/MMins.Fs;
MCSEM.K2O.Py=0.0*MOxide.K2O/MMins.Py;
MCSEM.K2O.Gs=0.0*MOxide.K2O/MMins.Gs;
MCSEM.K2O.Sp=0.0*MOxide.K2O/MMins.Sp;
MCSEM.K2O.Us=0.0*MOxide.K2O/MMins.Us;
MCSEM.K2O.An=0.0*MOxide.K2O/MMins.An;
MCSEM.K2O.Ab=0.0*MOxide.K2O/MMins.Ab;

%P2O5 mass fraction in all end-members
MCSEM.P2O5.Fo=0.0*MOxide.P2O5/MMins.Fo;
MCSEM.P2O5.Fa=0.0*MOxide.P2O5/MMins.Fa;
MCSEM.P2O5.Di=0.0*MOxide.P2O5/MMins.Di;
MCSEM.P2O5.Hd=0.0*MOxide.P2O5/MMins.Hd;
MCSEM.P2O5.En=0.0*MOxide.P2O5/MMins.En;
MCSEM.P2O5.Fs=0.0*MOxide.P2O5/MMins.Fs;
MCSEM.P2O5.Py=0.0*MOxide.P2O5/MMins.Py;
MCSEM.P2O5.Gs=0.0*MOxide.P2O5/MMins.Gs;
MCSEM.P2O5.Sp=0.0*MOxide.P2O5/MMins.Sp;
MCSEM.P2O5.Us=0.0*MOxide.P2O5/MMins.Us;
MCSEM.P2O5.An=0.0*MOxide.P2O5/MMins.An;
MCSEM.P2O5.Ab=0.0*MOxide.P2O5/MMins.Ab;

%H2O mass fraction in all end-members
MCSEM.H2O.Fo=0.0*MOxide.H2O/MMins.Fo;
MCSEM.H2O.Fa=0.0*MOxide.H2O/MMins.Fa;
MCSEM.H2O.Di=0.0*MOxide.H2O/MMins.Di;
MCSEM.H2O.Hd=0.0*MOxide.H2O/MMins.Hd;
MCSEM.H2O.En=0.0*MOxide.H2O/MMins.En;
MCSEM.H2O.Fs=0.0*MOxide.H2O/MMins.Fs;
MCSEM.H2O.Py=0.0*MOxide.H2O/MMins.Py;
MCSEM.H2O.Gs=0.0*MOxide.H2O/MMins.Gs;
MCSEM.H2O.Sp=0.0*MOxide.H2O/MMins.Sp;
MCSEM.H2O.Us=0.0*MOxide.H2O/MMins.Us;
MCSEM.H2O.An=0.0*MOxide.H2O/MMins.An;
MCSEM.H2O.Ab=0.0*MOxide.H2O/MMins.Ab;