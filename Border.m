function TFCout=Border(TFCin)
%This function is used to update and fill all borders of properties of system.
%The reason why I used this function to update properties is 
%there are many field properties that should be updated, 
%I do not want type them mannually.
%Created on 2020-7-3

global NIX
global NIY
    %TFCin(1,2:NIX+1)=TFCin(2,2:NIX+1);%top cool boundary solid stays the same
    TFCin(NIY+2,2:NIX+1)=TFCin(NIY+1,2:NIX+1);%bottom insulated boundary
    TFCin(1:NIY+2,1)=TFCin(1:NIY+2,2);%left insulated boundary
    TFCin(1:NIY+2,NIX+2)=TFCin(1:NIY+2,NIX+1);%right insulated boundary

    TFCout=TFCin;

end

