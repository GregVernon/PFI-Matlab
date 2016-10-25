function OutputData(cCOORD,pCOORD,RES,t,dCycle)

nVars = length(RES.VelProfile);
for n = 1:nVars
%     fileName = RES.VelProfile{n}.fileName;
    RES.VelProfile{n}.time = t;
    writeData = RES.VelProfile{n};
    dSize = size(RES.VelProfile{n}.DATA);
    if dCycle == 1
        RES.VelProfile{n}.file.DATA(1:dSize(1),1:dSize(2),1:dSize(3)) = RES.VelProfile{n}.DATA;
%         save(fileName,'-struct','writeData','-v7.3','-mat');
    else
        [nRows,nCols,nPages] = size(RES.VelProfile{n}.file,'DATA');
        ROWS = nRows+1 : nRows+dSize(1);
        COLS = nCols+1 : nCols+dSize(2);
        PAGES = nPages+1 : nPages+dSize(3);
        RES.VelProfile{n}.file.DATA(1:dSize(1),1:dSize(2),PAGES) = RES.VelProfile{n}.DATA;
%         save(fileName,'-struct','writeData','-v7.3','-mat');
    end
end