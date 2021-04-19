function [FRCPST,FKT,FR] = INTSP(FRTM,CPLTM,RLTM)
%1.To integrate solid properties, FS*RS*CpS, FS*RS*MCS, FS*RS*TCS, FS*KS, FS*RS; since
%KS, CpS, RS are functions of temperature and composition which is not homogeneous in solidified
%portion of finite volume.
%To sum over all minerals to give volume fraction-based effective
%properties. 
%2.To calculate current step RS*CpS, RL*CpL;
%Created 2019-10-17
%Modified for MAGTC3 2020-6-30

%---------- INPUT: old quantity -------
%FRCPSTM = FS*RS*CPS [NIY+2,NIX+2]
%FRCSMjTM = FS*RS*MCS [NIY+2,NIX+2]*major elements
%FRCSTeTM = FS*RS*TCS [NIY+2,NIX+2]*trace elements
%FKTM = FS*KS [NIY+2,NIX+2]
%FRTM = FS*RS [NIY+2,NIX+2]
%CPSTM [NIY+2,NIX+2]
%CPLTM [NIY+2,NIX+2]
%MCSTM [NIY+2,NIX+2]
%TCSTM [NIY+2,NIX+2]
%KSTM [NIY+2,NIX+2]
%RSTM [NIY+2,NIX+2]
%RLTM [NIY+2,NIX+2]
%dFSTM [NIY+2,NIX+2]

%---------- OUTPUT: new quantity -------
%FRCPS: FS*RS*CpS [NIY+2,NIX+2]
%FRCSMj: FS*RS*MCS [NIY+2,NIX+2]*major elements
%FRCSTe: FS*RS*TCS [NIY+2,NIX+2]*trace elements
%FK: FS*KS [NIY+2,NIX+2]
%FR: FS*RS [NIY+2,NIX+2]
%RSCPS: RS*CpS [NIY+2,NIX+2]
%RLCPL: RL*CpL [NIY+2,NIX+2]

global NIX
global NIY
global FRCSMj
global FRCSTe
global Majors
global Minors
global MCS
global TCS
global Mins
global dFS
global RS
global FRCPS
global CPS
global FK
global KS

FRCPST=zeros(NIY+2,NIX+2);%total solid FS*RS*CPS
FKT=zeros(NIY+2,NIX+2);%total solid FS*KS
FR=zeros(NIY+2,NIX+2);%FS*RS

%New quantity=Old quantity+delta(quantity)
%In fact, only one sentence is needed (matrix addition); but we unfold these formulae for
%better understanding as following
for i=1:NIX+2
    for j=1:NIY+2
        %------------------------- FORTRAN --------------------------------
        %FS*RS*CPS for solid
        FRCPS.OL(j,i)=FRCPS.OL(j,i)+dFS.OL(j,i)*RS.OL(j,i)*CPS.OL(j,i);
        FRCPS.OPX(j,i)=FRCPS.OPX(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*CPS.OPX(j,i);
        FRCPS.CPX(j,i)=FRCPS.CPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*CPS.CPX(j,i);
        FRCPS.PL(j,i)=FRCPS.PL(j,i)+dFS.PL(j,i)*RS.PL(j,i)*CPS.PL(j,i);
        FRCPS.ILM(j,i)=FRCPS.ILM(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*CPS.ILM(j,i);
        
        %FS*KS for solid
        FK.OL(j,i)=FK.OL(j,i)+dFS.OL(j,i)*KS.OL(j,i);
        FK.OPX(j,i)=FK.OPX(j,i)+dFS.OPX(j,i)*KS.OPX(j,i);
        FK.CPX(j,i)=FK.CPX(j,i)+dFS.CPX(j,i)*KS.CPX(j,i);
        FK.PL(j,i)=FK.PL(j,i)+dFS.PL(j,i)*KS.PL(j,i);
        FK.ILM(j,i)=FK.ILM(j,i)+dFS.ILM(j,i)*KS.ILM(j,i);
        
%         %total FS*RS*MCS for major elements in solid (This will be updated in MAGTFC.m)
%         FRCSMj.SiO2(j,i)=FRCSMj.SiO2(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.SiO2.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.SiO2.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.SiO2.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.SiO2.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.SiO2.ILM(j,i);
%         FRCSMj.TiO2(j,i)=FRCSMj.TiO2(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.TiO2.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.TiO2.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.TiO2.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.TiO2.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.TiO2.ILM(j,i);
%         FRCSMj.Al2O3(j,i)=FRCSMj.Al2O3(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.Al2O3.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.Al2O3.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.Al2O3.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.Al2O3.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.Al2O3.ILM(j,i);
%         FRCSMj.FeO(j,i)=FRCSMj.FeO(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.FeO.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.FeO.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.FeO.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.FeO.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.FeO.ILM(j,i);
%         FRCSMj.Fe2O3(j,i)=FRCSMj.Fe2O3(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.Fe2O3.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.Fe2O3.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.Fe2O3.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.Fe2O3.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.Fe2O3.ILM(j,i);
%         FRCSMj.MnO(j,i)=FRCSMj.MnO(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.MnO.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.MnO.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.MnO.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.MnO.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.MnO.ILM(j,i);
%         FRCSMj.MgO(j,i)=FRCSMj.MgO(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.MgO.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.MgO.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.MgO.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.MgO.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.MgO.ILM(j,i);
%         FRCSMj.CaO(j,i)=FRCSMj.CaO(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.CaO.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.CaO.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.CaO.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.CaO.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.CaO.ILM(j,i);
%         FRCSMj.Na2O(j,i)=FRCSMj.Na2O(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.Na2O.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.Na2O.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.Na2O.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.Na2O.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.Na2O.ILM(j,i);
%         FRCSMj.K2O(j,i)=FRCSMj.K2O(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.K2O.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.K2O.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.K2O.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.K2O.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.K2O.ILM(j,i);
%         FRCSMj.P2O5(j,i)=FRCSMj.P2O5(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.P2O5.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.P2O5.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.P2O5.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.P2O5.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.P2O5.ILM(j,i);
%         FRCSMj.H2O(j,i)=FRCSMj.H2O(j,i)+dFS.OL(j,i)*RS.OL(j,i)*MCS.H2O.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*MCS.H2O.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*MCS.H2O.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*MCS.H2O.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*MCS.H2O.ILM(j,i);
%         
%         %total FS*RS*TCS for minor elements in solid (This will be updated in MAGTFC.m)
%         FRCSTe.Sm(j,i)=FRCSTe.Sm(j,i)+dFS.OL(j,i)*RS.OL(j,i)*TCS.Sm.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*TCS.Sm.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*TCS.Sm.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*TCS.Sm.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*TCS.Sm.ILM(j,i);
%         FRCSTe.Nd(j,i)=FRCSTe.Nd(j,i)+dFS.OL(j,i)*RS.OL(j,i)*TCS.Nd.OL(j,i)+dFS.OPX(j,i)*RS.OPX(j,i)*TCS.Nd.OPX(j,i)+dFS.CPX(j,i)*RS.CPX(j,i)*TCS.Nd.CPX(j,i)+dFS.PL(j,i)*RS.PL(j,i)*TCS.Nd.PL(j,i)+dFS.ILM(j,i)*RS.ILM(j,i)*TCS.Nd.ILM(j,i);
        
        %FS*RS total density of solid
        FR(j,i)=FRTM(j,i)+RS.OL(j,i)*dFS.OL(j,i)+RS.OPX(j,i)*dFS.OPX(j,i)+RS.CPX(j,i)*dFS.CPX(j,i)+RS.PL(j,i)*dFS.PL(j,i)+RS.ILM(j,i)*dFS.ILM(j,i);
        
        %------------------------- MATLAB --------------------------------
%         %FS*RS*CPS for solid
%         for k=1:length(Mins)
%             cmd=['FRCPS.',Mins{k},'(j,i)=','FRCPS.',Mins{k},'(j,i)+dFS(j,i).',Mins{k},'(j,i)*RS.',Mins{k},'(j,i)*CPS.',Mins{k},'(j,i);'];
%             eval(cmd);
%         end
%         
%         %FS*KS for solid
%         for k=1:length(Mins)
%             cmd=['FK.',Mins{k},'(j,i)=FK.',Mins{k},'(j,i)+dFS.',Mins{k},'(j,i)*KS.',Mins{k},'(j,i);'];
%             eval(cmd);
%         end
%         
%         %FS*RS*MCS of major elements in solid
%         for k=1:length(Majors)
%             Chems=0.0;
%             for m=1:length(Mins)
%                 cmd=['dFS.',Mins{m},'(j,i)*RS.',Mins{m},'(j,i)*MCS.',Majors{k},'.',Mins{m},'(j,i)'];
%                 Chems=Chems+eval(cmd);
%             end
%             cmd=['FRCSMj.',Majors{k},'(j,i)=','FRCSMj.',Majors{k},'(j,i)+Chems'];
%             eval(cmd);
%         end
%         
%         %FS*RS*TCS of trace elements in solid
%         for k=1:length(Minors)
%             Chems=0.0;
%             for m=1:length(Mins)
%                 cmd=['dFS.',Mins{m},'(j,i)*RS.',Mins{m},'(j,i)*TCS.',Minors{k},'.',Mins{m},'(j,i)'];
%                 Chems=Chems+eval(cmd);
%             end
%             cmd=['FRCSTe.',Minors{k},'(j,i)=','FRCSTe.',Minors{k},'(j,i)+Chems'];
%             eval(cmd);
%         end
%         
%         %FS*RS total density of solid
%         Chems=0.0;
%         for k=1;length(Mins)
%             cmd=['RS.',Mins{k},'(j,i)*dFS.',Mins{k},'(j,i)'];
%             Chems=Chems+eval([cmd]);
%         end
%         FRT(j,i)=FRTM(j,i)+Chems;
        
        FRCPST(j,i)=FRCPS.OL(j,i)+FRCPS.OPX(j,i)+FRCPS.CPX(j,i)+FRCPS.PL(j,i)+FRCPS.ILM(j,i);
        FKT(j,i)=FK.OL(j,i)+FK.OPX(j,i)+FK.CPX(j,i)+FK.PL(j,i)+FK.ILM(j,i);
        
    end
end

end

