function OutputData(cCOORD,pCOORD,RES,t,dCycle)

nVars = length(RES.VelProfile);
for n = 1:nVars
    fileName = RES.VelProfile{n}.fileName;
    RES.VelProfile{n}.time = t;
    writeData = RES.VelProfile{n};
    if dCycle == 1
        save(fileName,'-struct','writeData','-v7.3','-mat');
    else
        save(fileName,'-struct','writeData','-v7.3','-mat');
    end
end