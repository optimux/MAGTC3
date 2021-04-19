function MUOE=VisMix(FSOTTM,MUOTM)
%This function calculates effective or bulk dynamic viscosity of mixture
%Created 2020-7-24

global FSCR2
global FSCR1
global NIX
global NIY

[r,c]=size(MUOTM);

%======= INPUTS =======
%FSOTTM: Total old FS [NIY+2,NIX+2]
%MUOTM: old pure liquid viscosity [NIY+2,NIX+2]

%======= OUTPUTS ======
%MUOE: effective or bulk dynamic viscosity of mixture [NIY+2,NIX+2]

%% ================= Einstein(1906) ====================
%Solid volume fraction should be small
%MUOE=(1.0+2.5*FSOTTM).*MUOTM;

%% ================== Roscoe(1952) =====================
%Set phi_m=FSCR2
%MUOE=MUOTM./(1.0-FSOTTM/FSCR2).^2.5;

%% ================ Solomatov(2007) ====================

%% =============== Costa(2005) Li2Si2O5 ================
%Eq.(9)
%% =============== Costa(2005) Pyrope ==================
%Eq.(9)
%% =============== Lesher 2015 ==================
%0<=FS<=FSCR1, free-moving crystals bulk viscosity [Pa.sec]
MUOE=zeros(r,c);
for i=1:c
    for j=1:r
        if(FSOTTM(j,i)<=FSCR1)
            MUOE(j,i)=MUOTM(j,i)*((1.0-FSOTTM(j,i))/(1.0-FSOTTM(j,i)/FSCR2))^((3.0*FSCR2)/(1.0-FSCR2));
        else
            MUOE(j,i)=MUOTM(j,i);
        end
    end
end

end

