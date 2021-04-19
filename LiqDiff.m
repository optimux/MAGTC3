function LiqDiff(TTM,PATM)
%This function is used to calculate diffusivity of oxides at given pressure and temperature, following
%equation 5-15 and Table S5.10 in Lesher&Spera2015
%Created on 2020-7-2

global NIX
global NIY
global DL0
global Ea
global Va
global DL

%TTM: temperature [K]
%PATM: absolute pressure [Pa]
%Ea: Activation energy [J/mol]
%Va: Activation volume [m^3/mol]
%R: universal gas constant 8.314 [J/mol/K]

%DL=DL0*exp(-(Ea+PVa)/RT)
for i=1:NIX+2
    for j=1:NIY+2
        DL.SiO2(j,i)=DL0.SiO2*exp(-(Ea.SiO2+PATM(j,i)*Va.SiO2)/(TTM(j,i)*8.314));
        DL.TiO2(j,i)=DL0.TiO2*exp(-(Ea.TiO2+PATM(j,i)*Va.TiO2)/(TTM(j,i)*8.314));
        DL.Al2O3(j,i)=DL0.Al2O3*exp(-(Ea.Al2O3+PATM(j,i)*Va.Al2O3)/(TTM(j,i)*8.314));
        DL.FeO(j,i)=DL0.FeO*exp(-(Ea.FeO+PATM(j,i)*Va.FeO)/(TTM(j,i)*8.314));
        DL.Fe2O3(j,i)=DL0.Fe2O3*exp(-(Ea.Fe2O3+PATM(j,i)*Va.Fe2O3)/(TTM(j,i)*8.314));
        DL.MnO(j,i)=DL0.MnO*exp(-(Ea.MnO+PATM(j,i)*Va.MnO)/(TTM(j,i)*8.314));
        DL.MgO(j,i)=DL0.MgO*exp(-(Ea.MgO+PATM(j,i)*Va.MgO)/(TTM(j,i)*8.314));
        DL.CaO(j,i)=DL0.CaO*exp(-(Ea.CaO+PATM(j,i)*Va.CaO)/(TTM(j,i)*8.314));
        DL.Na2O(j,i)=DL0.Na2O*exp(-(Ea.Na2O+PATM(j,i)*Va.Na2O)/(TTM(j,i)*8.314));
        DL.K2O(j,i)=DL0.K2O*exp(-(Ea.K2O+PATM(j,i)*Va.K2O)/(TTM(j,i)*8.314));
        DL.P2O5(j,i)=DL0.P2O5*exp(-(Ea.P2O5+PATM(j,i)*Va.P2O5)/(TTM(j,i)*8.314));
        DL.H2O(j,i)=DL0.H2O*exp(-(Ea.H2O+PATM(j,i)*Va.H2O)/(TTM(j,i)*8.314));
        DL.Sm(j,i)=DL0.Sm*exp(-(Ea.Sm+PATM(j,i)*Va.Sm)/(TTM(j,i)*8.314));
        DL.Nd(j,i)=DL0.Nd*exp(-(Ea.Nd+PATM(j,i)*Va.Nd)/(TTM(j,i)*8.314));
        
        % You can add more chemical diffusivity here!
        % DL.SiO2(j,i)=DL0.SiO2*exp(-(Ea.SiO2+PATM(j,i)*Va.SiO2)/(TTM(j,i)*8.314));
        % DL.SiO2(j,i)=DL0.SiO2*exp(-(Ea.SiO2+PATM(j,i)*Va.SiO2)/(TTM(j,i)*8.314));
    end
end

%NOTE: DL at ghost cells are set equal to inner layer since the
%medium surrouding this cavity is not given, hence DL is not given by
%T in ghost cells
DL.SiO2=Border(DL.SiO2);
DL.TiO2=Border(DL.TiO2);
DL.Al2O3=Border(DL.Al2O3);
DL.FeO=Border(DL.FeO);
DL.Fe2O3=Border(DL.Fe2O3);
DL.MnO=Border(DL.MnO);
DL.MgO=Border(DL.MgO);
DL.CaO=Border(DL.CaO);
DL.Na2O=Border(DL.Na2O);
DL.K2O=Border(DL.K2O);
DL.P2O5=Border(DL.P2O5);
DL.H2O=Border(DL.H2O);
DL.Sm=Border(DL.Sm);
DL.Nd=Border(DL.Nd);

end

