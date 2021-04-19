function KSTM=Conduct(TTM,PATM)
%This script is used to calculate total thermal conductivity of solid
%phases, i.e., lattice + radiation, up to 2000 K
%Created on 2020-6-28

%TTM: [K]
%PATM: [GPa]

PAG=PATM/10^9;%Absolute pressure in GPa
KSTM=struct('OL',zeros(1,1),...%olivine thermal conductivity [W/m/K]
    'OPX',zeros(1,1),...%opx thermal conductivity [W/m/K]
    'CPX',zeros(1,1),...%cpx thermal conductivity [W/m/K]
    'PL',zeros(1,1),...%pl thermal conductivity [W/m/K]
    'ILM',zeros(1,1));%ilmentite thermal conductivity [W/m/K]

%% ==================== Fayalite ========================
% a0=0.2386e-4;%thermal expansivity from Fei 1995, also Shutt&Lesher2006
% a1=1.153e-8;
% a2=-0.0518;
% KT=135.1;%Hofmeister 1999
% ga=1.45;%Hofmeister 1999
% dKTdP=4.0;%Hofmeister 1999
% K298=3.16;%Hofmeister 1999
% a=0.33;%Hofmeister 1999
% 
% Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0);%Eq.(10) in Hofmeister 1999
% 
% KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
% KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals, Hofmeister 1999
% KSTM.OL=KLA+KRA;%lattice+radiation

%% ==================== Forsterite ========================
a0=0.3407e-4;%thermal expansivity from Fei 1995, also Shutt&Lesher2006
a1=0.8674e-8;
a2=-0.7545;

K298=5.2;%Data from Hofmeister 1999 Table 1
a=0.45;%fitting parameter
dKTdP=4.0;%Hofmeister 1999
KT=127.9;%isothermal bulk modulus
ga=1.25;%Gruneisen parameter

Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0);%Eq.(10) in Hofmeister 1999
KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
KSTM.OL=KLA+KRA;%lattice+radiation

%% ==================== Ferrosilite ========================
% a0=0.308e-4;%thermal expansivity from Shutt&Lesher2006
% a1=0.978e-8;
% a2=-0.404;
% a3=1.52e-18;
% 
% RHO=densityfull(298,1.0e-4);%reference at 298 K, 1 bar for K298
% %Fs=Ferrosilite
% Fs=178.7-0.001378*298-355550.0/298^2-1496.3/sqrt(298);%[J/mol/K]
% Fs=2.0*Fs/0.2638614;%Fe2Si2O6: 263.8614 g/mol [J/kg/K]
% Dlat=0.5e-6-(0.01753-0.00010365*298+2.2451e-7*298^2-3.407e-11*298^3)/(Fs*RHO{17,2});%D_total = 0.5e-6 m^2/sec at 298 K from Fig.3 of Hofmeister2014
% K298=Dlat*Fs*RHO{17,2};%17: ortho-Fs
% a=0.33;%fitting parameter, Hofmeister1999
% dKTdP=9.0;%Hofmeister 1999
% KT=97.864;%isothermal bulk modulus
% %KT=98.8/(1+1.1*(0.308e-4+298*0.978e-8-0.404/298^2+1.52e-18*298^4)*298)
% ga=1.1;%Gruneisen parameter
% 
% 
% Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0);%Eq.(10) in Hofmeister 1999
% KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
% KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals, Hofmeister 1999
% KSTM.OPX=KLA+KRA;%lattice+radiation

%% ==================== Enstatite ========================
%[1]: Hofmeister1999
a0=0.271e-4;%thermal expansivity from Shutt&Lesher2006
a1=1.2e-8;
a2=-0.66;
a3=3.1e-18;
dKTdP=10.2;%Schutt&Lesher2006, also conforms Table 1 in Hofmeister 1999
KT=107;%Hofmeister1999
ga=0.96;%Hofmeister1999
a=0.33;

Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0)+0.2*a3*(TTM^5-298^5);%Caption of Table 3 in Schutt&Lesher2006
KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
KSTM.OPX=KLA+KRA;

%[2]: Hofmeister2012 Table 2
%KLA=1.0e-6/(-0.67915+0.0051741*TTM-3.6072e-6*TTM^2+1.0339e-9*TTM^3);%up to 1275 K
%Note[1]: partial(ln(k))/partial(P) ~0.04 per GPa
%Note[2]: no way to figure out ferrosilite


%% ==================== Hedenbergite ========================
% a0=0.232e-4;%thermal expansivity from Shutt&Lesher2006
% a1=1.88e-8;
% a2=0.0;
% 
% if(TTM<=1100.0)%max fit for T 1100 K
%     Dlat=1.0e-6/(-0.36581+0.0040009*TTM-15.671e-7*TTM^2);%[m^2/sec]
%     Hd=155.2+0.006285*TTM-923000.0/TTM^2-1020.0/sqrt(TTM);%[J/mol/K]
%     Hd=2.0*Hd/0.2480924;%CaFeSi2O6: 248.0924 g/mol [J/kg/K]
%     RHO=densityfull(TTM,PAG);%input K, GPa
%     KLA=Dlat*Hd*RHO{30,2};%30: clino-Hd
%     KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals like olvine, Hofmeister 1999
%     %Assuming Aegirine NaFe(III)Si2O6 as Hd end-member, its lattice diffusivity up to 1100 K
%     KSTM.CPX=KRA+KLA;
% else
%     KLA=1.7478;%KLA(1100K)=1.7478 [m^2/sec]
%     KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals like olvine, Hofmeister 1999
%     %Assuming Aegirine NaFe(III)Si2O6 as Hd end-member, its lattice diffusivity up to 1100 K
%     KSTM.CPX=KRA+KLA;
% end

%Note: pressure derivative partial(ln(k_lat))/partial(P) ~ 0.042-0.047 per GPa

%% ==================== Diopside ========================
%[1]: K(T) is replaced by poly fit of 1/D in Table 2 of Hofmeister 2008
a0=0.232e-4;%thermal expansivity from Shutt&Lesher2006
a1=1.88e-8;
a2=0.0;
dKTdP=4.8;%Schutt&Lesher 2006
KT=109.56;%derived from Schutt&Lesher 2006
ga=1.0;%Schutt&Lesher 2006


if(TTM<=1400.0)%Lattice saturation tempearture 1400 K
    D=1.0e-6/(-0.22631+0.0023854*TTM-6.9593e-7*TTM^2);%thermal diffusivity from Hofmeister 2008 Table 2 [m^2/sec]
    %NOTE: we assume this diopside is pure Di from Hofmerister 2008
    Di=157.3+0.0000205*TTM-1372950.0/TTM^2-1010.5/sqrt(TTM);%[J/mol/K]
    Di=2.0*Di/0.2165504;%CaMgSi2O6: 216.5504 g/mol [J/kg/K]
    RHO=densityfull(TTM,PAG);%input K, GPa
    KLA=Di*D*RHO{29,2};%29: clino-Di
    KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999.
    KSTM.CPX=KLA+KRA;
else
    Di=157.3+0.0000205*TTM-1372950.0/TTM^2-1010.5/sqrt(TTM);%[J/mol/K]
    Di=2.0*Di/0.2165504;%CaMgSi2O6: 216.5504 g/mol [J/kg/K]
    RHO=densityfull(TTM,PAG);%input K, GPa
    KLA=Di*5.7e-7*RHO{29,2};%29: clino-Di
    KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999.
    KSTM.CPX=KLA+KRA;
end

%% ==================== Albite ========================
% %[1]: Albite from Hoemeister 2009 Table 6, only k_lat
% if(TTM<=1400.0)%Tg=1260 K, melting point 1400 K
%     %Albite [010],[001],[100] mean value as D, data from Table 6 of Hofmeister 2009
%     D=1.0/(0.493+257.26/TTM+27922.0/TTM^2)+...
%         1.0/(0.31095+314.91/TTM)+...
%         1.0/(0.33964+105.38/TTM+30164.0/TTM^2);%mm^2/sec
%     D=1.0e-6/(D/3.0);%m^2/sec
%     %Note: bulk value is calculated following Eq.(6) in Branlund 2012
%     
%     %It's Analbite
%     Ab=309.74+0.015278*TTM-26160000.0/TTM^2+4109100000.0/TTM^3+8840.0/TTM;%[J/mol/K]
%     Ab=Ab/0.262223;%NaAlSi3O8: 262.223 g/mol [J/kg/K]
%     RHO=densityfull(TTM,PAG);%input K, GPa
%     KLA=D*Ab*RHO{46,2};%46: Albite
%     KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999.
%     KSTM.PL=KRA+KLA;
% else
%     D=0.5e-6;%Dsat,melt
%     Ab=309.74+0.015278*TTM-26160000.0/TTM^2+4109100000.0/TTM^3+8840.0/TTM;%[J/mol/K]
%     Ab=Ab/0.262223;%NaAlSi3O8: 262.223 g/mol [J/kg/K]
%     RHO=densityfull(TTM,PAG);%input K, GPa
%     KLA=D*Ab*RHO{46,2};%46: Albite
%     KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999.
%     KSTM.PL=KRA+KLA;
% end

%% ==================== Sanidine ========================
% %[1]: KAlSi3O8(sanidine) from Pertermann 2008, only k_lat
% %sanidine Tg=1370 K, melting point 1400 K
% 
% %OFD [100],[010],[001] mean value as D_bulk, data from Table 4 of Pertermann2008
% D_bulk=1.0/(0.4183+82.33/TTM+9356.8/TTM^2)+...
%     1.0/(0.699+99.88/TTM+25726.6/TTM^2)+...
%     1.0/(0.5492+110.4/TTM+15574.4/TTM^2);%mm^2/sec
% D_bulk=3.0/(D_bulk*1.0e-6);%m^2/sec
% %Note: bulk value is calculated following Eq.(6) in Branlund 2012
% 
% D_melt=0.52615-9.1562/TTM+9920.3/TTM^2;%mm^2/sec
% D_melt=D_melt*1.0e-6;%m^2/sec

%NOTE: no CPS of sanidine, so omitted

%% ==================== Anorthite ========================
b0=2.44e-5-3.1e-7*1.0+1.8e-9*1.0^2;%thermal expansivity from Korenaga&Keronaga2016, X_an=1.0
b1=9.0e-9-4.0e-11*1.0;

%NOTE: both Eq.(6) and FMA in Table 2 of Branlund 2012 are not applicable
%here, since data of FMA is not continuous and scattered. So we use fitting
%of Eq.(7) in Christopher J. Grose and Juan Carlos Afonso 2013.
%Thank god, I read this paper! Even I have no idea why I read it!
%In fact, Grose&Afonso2013 also use FMA in Table 2 of Branlund 2012! 2020-1-24
D=(0.36+0.4*exp(-(TTM-273.0)/300.0))*1.0e-6;%[m^2/sec]

b=b0+2.0*b1*(TTM-298.0);
RHO=densityfull(TTM,PAG);%input K, GPa
An=290.90+0.0276*TTM-34080000.0/TTM^2+5218000000.0/TTM^3+29625.0/TTM;%[J/mol/K]
An=An/0.2782073;%CaAl2Si2O8: 278.2073 g/mol [J/kg/K]
KLA=D*An*RHO{47,2};%47: Anorthite

KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
KSTM.PL=KLA+KRA;

%% ==================== Grossular ========================
% %Hofmeister2006 Table 3
%         if(TTM<=1360.0)%fitting range
%             D=1.0e-6/(-0.10641+0.0017914*TTM-6.4912e-7*TTM^2);%[m^2/sec]
%
%             %Grossular
%             Gs=542.6+0.01294*TTM-3186000.0/TTM^2+277700000.0/TTM^3-56020.0/TTM;%[J/mol/K]
%             Gs=Gs/0.450446376;%Ca3Al2Si3O12: 450.446376 g/mol [J/kg/K]
%             RHO=densityfull(TTM,PAG);%input K, GPa
%             KLA=D*Gs*RHO{37,2};%37: Grossular
%             KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
%             KSTM.GRT=KLA+KRA;
%         end
%
% %% ==================== Andradite ========================
% %Hofmeister2006 Table 3
%         if(TTM<=573.0)%fitting range
%             D=1.0e-6/(-0.26309+0.0024971*TTM-1.4377e-6*TTM^2);%[m^2/sec]
%
%             %Andradite set as Gs
%             Gs=542.6+0.01294*TTM-3186000.0/TTM^2+277700000.0/TTM^3-56020.0/TTM;%[J/mol/K]
%             Gs=Gs/0.450446376;%Ca3Al2Si3O12: 450.446376 g/mol [J/kg/K]
%             RHO=densityfull(TTM,PAG);%input K, GPa
%             KLA=D*Gs*RHO{36,2};%36: Andradite
%             KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals, Hofmeister 1999
%             KSTM.GRT=KLA+KRA;
%         end

%% ==================== Pyrope ========================
% %Thermal properties from Schutt2006
% a0=0.2311e-4;
% a1=0.5956e-8;
% a2=-0.4538;
% dKTdP=1.5;%Schutt2006
% %dKTdP=3.8;%Hofmeister1999
% %KT=171.5;%Hofmeister1999
% KT=170.8;%Derived from Schutt2006
% %KT=172.1/(1+1.29*(0.2311e-4+298*0.5956e-8-0.4538/298^2)*298)
% %ga=1.43;%Hofmeister1999
% ga=1.29;%Schutt2006
% K298=3.2;%Hofmeister1999 == Horai1971
% a=0.33;
%
%         Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0);%Eq.(10) in Hofmeister 1999
%         KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
%         KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
%         KSTM.GRT=KLA+KRA;%lattice+radiation
%
%

% %% ==================== Magnetite ========================
% if(TTM<=500.0)%fitting range
%     D=1.0e-6/(-0.56653+0.0046202*TTM-3.4391e-6*TTM^2);%[m^2/sec]
%     
%     Mag=-3.558e3+3.3473e2*sqrt(TTM)-9.309*TTM+2.5388e-3*TTM^2+1.4273e5/TTM;%[J/mol/K] from Robie 1978
%     Mag=Mag/0.2315386;%Fe3O4: 231.5386 g/mol
%     RHO=densityfull(TTM,PAG);%input K, GPa
%     KLA=D*Mag*RHO{43,2};%43: Magnetite
%     KRA=0.01753-0.00010365*TTM+2.2451e-7*TTM^2-3.407e-11*TTM^3;%radiation only, for Ferrous minerals, Hofmeister 1999
%     KSTM.SP=KLA+KRA;
% end
% 
% %% ==================== Spinel ========================
% %[1]: Hofmeister1999
% %Thermal properties from Schutt2006
% a0=0.187e-4;
% a1=0.975e-8;
% a2=-0.365;
% dKTdP=5.05;%Schutt2006
% %dKTdP=4.9;%Hofmeister1999
% %KT=195.2;%Hofmeister1999
% KT=196.8;%Derived from Schutt2006
% %KT= 198/(1+1.1*(0.187e-4+298*0.975e-8-0.265/298^2)*298)
% %ga=1.4;%Hofmeister1999
% ga=1.1;%Schutt2006
% K298=9.5;%Hofmeister1999 == Horai1971 == Schon2015
% a=0.33;
%
%         Tintg=a0*(TTM-298)+0.5*a1*(TTM^2-298^2)+a2*(-1.0/TTM+1.0/298.0);%Eq.(10) in Hofmeister 1999
%         KLA=K298*(298.0/TTM)^a*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*PAG/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
%         KRA=8.5e-11*TTM^3;%radiation only, for Fe-free silicates and oxides, Hofmeister 1999
%         KSTM.SP=KLA+KRA;%lattice+radiation

%% ==================== Ilmenite ========================
%Heimo2019, Sample T1, total thermal conductivity
if(TTM<=1360.0)
    KSTM.ILM= -2.7e-12*TTM^4+8.3e-9*TTM^3-7.4e-6*TTM^2+1.3e-3*TTM+2.3;%[W/m/K]
end
