function MCSTM=SYSEOS(EndMems)
%This function is Equation of State (EOS) of System, to calculate property of new time step:
%how much will crystallize of each mineral or their relative ratios, chemical composition of
%crystallized minerals. In thoery, this function can be replaced by alphamenlts package, but
%alphamelts is not that much reliable for non-calibrated composition. While for Di-An system,
%the amount of each mineral is rather well characterized, so we mainly deal with this case
%while other minerals are assumed known.
%If composition, say MORB or gabbro, is well calibrated, we shall use alphamelts instead, for
%example, when dealing with Vanadic Titanomagnetite deposit in Panzhihua, Sichuan.

%created on 2020-7-5

global MCSEM
global MCSN
global NIX
global NIY

EndMems='Di';
%NOTE: temporarily, we only have Di-An system, so MCSN is simplified here.
for i=1:NIX+2
    for j=1:NIY+2
        MCSN.SiO2.CPX(j,i)=MCSEM.SiO2.Di;
        MCSN.SiO2.PL(j,i)=MCSEM.SiO2.An;
        
        MCSN.MgO.CPX(j,i)=MCSEM.MgO.Di;
        MCSN.MgO.PL(j,i)=MCSEM.MgO.An;
        
        MCSN.CaO.CPX(j,i)=MCSEM.CaO.Di;
        MCSN.CaO.PL(j,i)=MCSEM.CaO.An;
        
        MCSN.Al2O3.CPX(j,i)=MCSEM.Al2O3.Di;
        MCSN.Al2O3.PL(j,i)=MCSEM.Al2O3.An;
    end
end
MCSTM=20200705;
end
