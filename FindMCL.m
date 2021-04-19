function FindMCL(X)
%This function finds interpolation of elements concentration
%Created on 2020-4-12

global MCL
global melt
global Majors
global NIX
global NIY

for i=1:NIX
    for j=1:NIY
        [dif,No]=min(abs(melt.MgO-X(j,i)));%X is MgO concentration interested
        Mg=[melt.MgO(No-1) melt.MgO(No) melt.MgO(No+1)];%nearest MgO in melt
        
        for k=1:length(Majors)
            cmd=['MC=[melt.',Majors{k},'(No-1) melt.',Majors{k},'(No) melt.',Majors{k},'(No+1)];'];%nearest major oxides concentration
            eval(cmd);
            cmd=['MCL.',Majors{k},'(j,i)=interp1(Mg,MC,X(j,i),''spline'');'];%cubic spline interpolation of property X
            eval(cmd);
        end
    end
end

end

