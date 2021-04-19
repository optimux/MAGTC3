function FindNSE(FSO,FSN,step)
%Original: [NFSS,NFSE] = FindNSE(dFSO,dFSN,step)
%To find step numbers which show solidification starts and ends,
%repectively, for step>=2
%Created 2019-10-20
%Modified for MAGTC3 2020-6-30

%NFSS: record step at which minerals (OL, OPX, CPX, PL, ILM) begins to form
%NFSS=struct('OL',ones(NIY+2,NIX+2,'int16'),...
%     'OPX',ones(NIY+2,NIX+2,'int16'),...
%     'CPX',ones(NIY+2,NIX+2,'int16'),...
%     'PL',ones(NIY+2,NIX+2,'int16'),...
%     'ILM',ones(NIY+2,NIX+2,'int16'));
%NFSE: record step at which no more minerals (OL, OPX, CPX, PL, ILM) to form
% NFSE=struct('OL',ones(NIY+2,NIX+2,'int16'),...
%     'OPX',ones(NIY+2,NIX+2,'int16'),...
%     'CPX',ones(NIY+2,NIX+2,'int16'),...
%     'PL',ones(NIY+2,NIX+2,'int16'),...
%     'ILM',ones(NIY+2,NIX+2,'int16'));

%NOTE: these two variables above are used to cope with inhomogeneous Cu
%distribution in solid

%dFSO: old dFS
%dFSN: new dFS
%FSO: old FS
%FSN: new FS

global NIX
global NIY
global NFSS
global NFSE

errFS=0.0001;

%======================= VERY IMPORTANT NOTE ==============================
%Sometimes, remelting may occur, thus dFS<0.0; and also somtimes during
%solidification process, there is no crystallization for some steps. So,
%the judgement below based on dFS is no better than the one based on FS!
%======================= VERY IMPORTANT NOTE ==============================

%% ================== OLD ALGORITHM dFS ==========================
% for i=1:NIX
%     for j=1:NIY
%         %Situation One: no solidification yet
%         if((dFSO(j,i)<=0.0)&&(dFSN(j,i)<=0.0))
%             NFSS(j,i)=step;%solidification start step number should change until real solid shows up
%             NFSE(j,i)=step;%solidification end step number should >=NFSS
%         end
%         %Situation Two: solidification just begins right now
%         if((dFSO(j,i)<=0.0)&&(dFSN(j,i)>0.0))
%             NFSS(j,i)=step;%solidification start step number should change until real solid shows up
%             NFSE(j,i)=step;%solidification end step number should >=NFSS
%         end
%         %Situation Three: solidification keeps working
%         if((dFSO(j,i)>0.0)&&(dFSN(j,i)>0.0))
%             %NFSS(j,i)=step;%solidification start step number has been recorded
%             NFSE(j,i)=step;%solidification end step number should >=now
%         end
%         %Situation Four: solidification ends
%         if((dFSO(j,i)>0.0)&&(dFSN(j,i)<=0.0))
%             %NFSS(j,i)=step;%solidification start step number has been recorded
%             NFSE(j,i)=step;%solidification end step number should >=now
%         end
%
%     end
% end

%% ===================== NEW ALGORITHM FS ========================
for i=1:NIX+2
    for j=1:NIY+2
        %Situation One: no solidification yet
        if((abs(FSO.OL(j,i))<=errFS)&&(abs(FSN.OL(j,i))<=errFS))
            NFSS.OL(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.OL(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.OPX(j,i))<=errFS)&&(abs(FSN.OPX(j,i))<=errFS))
            NFSS.OPX(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.OPX(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.CPX(j,i))<=errFS)&&(abs(FSN.CPX(j,i))<=errFS))
            NFSS.CPX(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.CPX(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.PL(j,i))<=errFS)&&(abs(FSN.PL(j,i))<=errFS))
            NFSS.PL(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.PL(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.ILM(j,i))<=errFS)&&(abs(FSN.ILM(j,i))<=errFS))
            NFSS.ILM(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.ILM(j,i)=step;%solidification end step number should >=NFSS
        end

        
        %Situation Two: solidification just begins right now
        if((abs(FSO.OL(j,i))<=errFS)&&(abs(FSN.OL(j,i))>errFS)&&(abs(FSN.OL(j,i)-1.0)>=errFS))
            NFSS.OL(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.OL(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.OPX(j,i))<=errFS)&&(abs(FSN.OPX(j,i))>errFS)&&(abs(FSN.OPX(j,i)-1.0)>=errFS))
            NFSS.OPX(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.OPX(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.CPX(j,i))<=errFS)&&(abs(FSN.CPX(j,i))>errFS)&&(abs(FSN.CPX(j,i)-1.0)>=errFS))
            NFSS.CPX(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.CPX(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.PL(j,i))<=errFS)&&(abs(FSN.PL(j,i))>errFS)&&(abs(FSN.PL(j,i)-1.0)>=errFS))
            NFSS.PL(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.PL(j,i)=step;%solidification end step number should >=NFSS
        end
        if((abs(FSO.ILM(j,i))<=errFS)&&(abs(FSN.ILM(j,i))>errFS)&&(abs(FSN.ILM(j,i)-1.0)>=errFS))
            NFSS.ILM(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE.ILM(j,i)=step;%solidification end step number should >=NFSS
        end

        
        %Situation Three: solidification keeps working
        if((abs(FSO.OL(j,i))>errFS)&&(abs(FSO.OL(j,i)-1.0)>=errFS)&&(abs(FSN.OL(j,i))>errFS)&&(abs(FSN.OL(j,i)-1.0)>=errFS))
            %NFSS.OL(j,i)=step;%solidification start step number has been recorded
            NFSE.OL(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.OPX(j,i))>errFS)&&(abs(FSO.OPX(j,i)-1.0)>=errFS)&&(abs(FSN.OPX(j,i))>errFS)&&(abs(FSN.OPX(j,i)-1.0)>=errFS))
            %NFSS.OPX(j,i)=step;%solidification start step number has been recorded
            NFSE.OPX(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.CPX(j,i))>errFS)&&(abs(FSO.CPX(j,i)-1.0)>=errFS)&&(abs(FSN.CPX(j,i))>errFS)&&(abs(FSN.CPX(j,i)-1.0)>=errFS))
            %NFSS.CPX(j,i)=step;%solidification start step number has been recorded
            NFSE.CPX(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.PL(j,i))>errFS)&&(abs(FSO.PL(j,i)-1.0)>=errFS)&&(abs(FSN.PL(j,i))>errFS)&&(abs(FSN.PL(j,i)-1.0)>=errFS))
            %NFSS.PL(j,i)=step;%solidification start step number has been recorded
            NFSE.PL(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.ILM(j,i))>errFS)&&(abs(FSO.ILM(j,i)-1.0)>=errFS)&&(abs(FSN.ILM(j,i))>errFS)&&(abs(FSN.ILM(j,i)-1.0)>=errFS))
            %NFSS.ILM(j,i)=step;%solidification start step number has been recorded
            NFSE.ILM(j,i)=step;%solidification end step number should >=now
        end

        
        %Situation Four: solidification ends
        if((abs(FSO.OL(j,i))>errFS)&&(abs(FSO.OL(j,i)-1.0)>=errFS)&&(abs(FSN.OL(j,i)-1.0)<errFS))
            %NFSS.OL(j,i)=step;%solidification start step number has been recorded
            NFSE.OL(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.OPX(j,i))>errFS)&&(abs(FSO.OPX(j,i)-1.0)>=errFS)&&(abs(FSN.OPX(j,i)-1.0)<errFS))
            %NFSS.OPX(j,i)=step;%solidification start step number has been recorded
            NFSE.OPX(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.CPX(j,i))>errFS)&&(abs(FSO.CPX(j,i)-1.0)>=errFS)&&(abs(FSN.CPX(j,i)-1.0)<errFS))
            %NFSS.CPX(j,i)=step;%solidification start step number has been recorded
            NFSE.CPX(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.PL(j,i))>errFS)&&(abs(FSO.PL(j,i)-1.0)>=errFS)&&(abs(FSN.PL(j,i)-1.0)<errFS))
            %NFSS.PL(j,i)=step;%solidification start step number has been recorded
            NFSE.PL(j,i)=step;%solidification end step number should >=now
        end
        if((abs(FSO.ILM(j,i))>errFS)&&(abs(FSO.ILM(j,i)-1.0)>=errFS)&&(abs(FSN.ILM(j,i)-1.0)<errFS))
            %NFSS.ILM(j,i)=step;%solidification start step number has been recorded
            NFSE.ILM(j,i)=step;%solidification end step number should >=now
        end
        
    end
end


end

