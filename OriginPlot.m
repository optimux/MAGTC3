function OriginPlot(VXP,VYP,TP)
%This function is used to plot in OriginLab format
%Created on 2020-7-3

global NIX
global NIY
global x
global y

VXPP=0.0;%liquid x-axis velocity [m/sec]
VYPP=0.0;%liquid y-axis velocity [m/sec]
VPP=0.0;%liquid absolute magnitute velocity [m/sec]
DIRL=0.0;%liquid velocity direction assuming 0.0 is right horizon arrow [radians]


Frame=zeros(NIY*NIX,4);

for i=1:NIX
    for j=1:NIY
        %--------------------- LIQUID ------------------------
        VXPP=0.5*(VXP(j+1,i)+VXP(j+1,i+1));%x-axis liquid velocity [m/sec]
        VYPP=0.5*(VYP(j,i+1)+VYP(j+1,i+1));%y-axis liquid velocity [m/sec]
        if((abs(VYPP)<=1.0e-9)&&(abs(VXPP)>1.0e-9))%VYPP=0.0, VXPP!=0.0
            if(VXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VXPP)<=1.0e-9)&&(abs(VYPP)>1.0e-9))%VXPP=0.0, VYPP!=0.0
            if(VYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VXPP)<=1.0e-9)&&(abs(VYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VXPP>0.0)&&(VYPP>0.0))
        DIRL=atan(VYPP/VXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VXPP>0.0)&&(VYPP<0.0))
            DIRL=atan(VYPP/VXPP)+2.0*pi;
        end
        if((VXPP<0.0)&&(VYPP>0.0))
            DIRL=atan(VYPP/VXPP)+pi;
        end
        if((VXPP<0.0)&&(VYPP<0.0))
            DIRL=pi+atan(VYPP/VXPP);
        end
        
                VPP=sqrt(VXPP^2+VYPP^2);%absolute magnitute of liquid velocity [m/sec]

        Frame(NIY*(i-1)+j,1)=x(i+1);
        Frame(NIY*(i-1)+j,2)=y(j+1);
        Frame(NIY*(i-1)+j,3)=DIRL;
        Frame(NIY*(i-1)+j,4)=VPP;
        
    end
end
xlswrite('VSV.xlsx',Frame(:,1:4),'LiqMov');
xlswrite('VSV.xlsx',TP,'T');
end

