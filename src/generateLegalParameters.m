function pGraph = generateLegalParameters
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);

INFORMATION.Children = [];
INFORMATION.Variables = struct('title',blockTitle,'description',description);
INFORMATION.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
nx = struct('name','nx','class',{'double'},'isRequired',true);
ny = struct('name','ny','class',{'double'},'isRequired',true);
xmin = struct('name','xmin','class',{'double'},'isRequired',true);
xmax = struct('name','xmax','class',{'double'},'isRequired',true);
ymin = struct('name','ymin','class',{'double'},'isRequired',true);
ymax = struct('name','ymax','class',{'double'},'isRequired',true);
INITIAL_GEOMETRY_DEFINITION.Children = [];
INITIAL_GEOMETRY_DEFINITION.Variables = struct('title',blockTitle,'description',description,...
                                               'nx',nx,'ny',ny,...
                                               'xmin',xmin,'xmax',xmax,...
                                               'ymin',ymin,'ymax',ymax);
INITIAL_GEOMETRY_DEFINITION.isRequired = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
relTol = struct('name','relative tolerance','class',{'double'},'isRequired',false,'default',1e-5);
solver = struct('name','solver','class','string','isRequired',false,'default','Jacobi','validValues',{{'Jacobi';'PCG';'GMRES';'Direct';'Linfactor';'Inverse'}});
precond = struct('name','preconditioner','class','string','isRequired',false,'default',[],'validValues',{{'ilu'}});
kspace = struct('name','k-space','class','double','isRequired',false,'default',20);
maxit = struct('name','maximum iterations','class','double','isRequired',false,'default',1e4);

LINEAR_SOLVE.Children = [];
LINEAR_SOLVE.Variables = struct('title',blockTitle,'description',description,...
                                'relTol',relTol,'solver',solver,...
                                'precond',precond,'kspace',kspace,...
                                'maxit',maxit);
LINEAR_SOLVE.isRequired = true;
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
stepTime = struct('name','step time','class','double','isRequired',true);
maxTime = struct('name','max time control','class','string','isRequired',false,'default','Reynolds');
scale = struct('name','max time scale','class','double','isRequired',false,'default',0.01);

TIME_CONTROL.Children = [];
TIME_CONTROL.Variables = struct('title',blockTitle,'description',description,...
                                'stepTime',stepTime,'maxTime',maxTime,...
                                'scale',scale);
TIME_CONTROL.isRequired = true;
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);

SOLUTION_CONTROL.Children = struct('LINEAR_SOLVE',LINEAR_SOLVE,'TIME_CONTROL',TIME_CONTROL);
SOLUTION_CONTROL.Variables = struct('title',blockTitle,'description',description);
SOLUTION_CONTROL.isRequired = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
func = struct('name','function type','class','string','isRequired',true,'validValues',{'constant';},'default',{'constant'});
value = struct('name','function value','class','double','isRequired',true);
scale = struct('name','scale factor','class','double','isRequired',false,'default',1);

REYNOLDS.Children = [];
REYNOLDS.Variables = struct('title',blockTitle,'description',description,...
                            'func',func,'value',value,'scale',scale);
REYNOLDS.isRequired = true;
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
func = struct('name','function type','class','string','isRequired',true,'validValues',{'constant';},'default',{'constant'});
value = struct('name','function value','class','double','isRequired',true);
scale = struct('name','scale factor','class','double','isRequired',false,'default',1);

CIRCULATION.Children = [];
CIRCULATION.Variables = struct('title',blockTitle,'description',description,...
                               'func',func,'value',value,'scale',scale);
CIRCULATION.isRequired = true;
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
FLOW_DEFINITION.Children = struct('REYNOLDS',REYNOLDS,'CIRCULATION',CIRCULATION);
FLOW_DEFINITION.Variables = struct('title',blockTitle,'description',description);
FLOW_DEFINITION.isRequired = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
outputName = struct('name','output filename','class','string','isRequired',true);
activeStep = struct('name','active step','class',{{'double';'string'}},'isRequired',false,'default','all');
outInit = struct('name','output initial step','class','logical','isRequired',false,'default',true);
outFin = struct('name','output final step','class','logical','isRequired',false,'default',true);
tIncr = struct('name','time increment','class','double','isRequired',true);
dBuffer = struct('name','buffer increment','class','double','isRequired',false,'default',100);
variable = struct('name','output variable','class','string','isRequired',true,'options',{{'OMEGA','PSI','pU','pV'}});
iX = struct('name','iX','class','double','isRequired',true);
iY = struct('name','jY','class','double','isRequired',true);

RESULT.Children = [];
RESULT.Variables = struct('title',blockTitle,'description',description,...
                          'outputName',outputName,'activeStep',activeStep,...
                          'outInit',outInit,'outFin',outFin,'tIncr',tIncr,...
                          'dBuffer',dBuffer,'variable',variable,...
                          'iX',iX,'iY',iY);
RESULT.isRequired = false;

%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
logName = struct('name','log filename','class','string','isRequired',false,'default',strcat('PFI_log_',datestr(now,'mm_dd_yyyy__HH:MM:SS')));
tIncr = struct('name','time increment','class','double','isRequired',false,'default',1e3);
initStep = struct('name','initial step','class','logical','isRequired',false,'default',true);
finStep = struct('name','final step','class','logical','isRequired',false,'default',true);
logData = struct('name','log data','class','cell','isRequired',false,'default',{{'Timestep';'Time';'Time Elapsed'}});

LOG_INFO.Children = [];
LOG_INFO.Variables = struct('title',blockTitle,'description',description,...
                            'logName',logName,'tIncr',tIncr,...
                            'initStep',initStep,'finStep',finStep,...
                            'logData',logData);
LOG_INFO.isRequired = false;

%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);
plotFlow = struct('name','plotFlow','class','logical','isRequired',false,'default',false);
plotState = struct('name','plotState','class','logical','isRequired',false,'default',false);
saveFlow = struct('name','saveFlow','class','logical','isRequired',false,'default',false);
saveState = struct('name','saveState','class','logical','isRequired',false,'default',false);

IN_SITU.Children = [];
IN_SITU.Variables = struct('title',blockTitle,'description',description,...
                           'plotFlow',plotFlow,'plotState',plotState,...
                           'saveFlow',saveFlow,'saveState',saveState);
IN_SITU.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
blockTitle = struct('name','title','class','string','isRequired',false);
description = struct('name','description','class','string','isRequired',false);

STEP.Children = struct('SOLUTION_CONTROL',SOLUTION_CONTROL,...
                       'FLOW_DEFINITION',FLOW_DEFINITION,...
                       'RESULT',RESULT,'LOG_INFO',LOG_INFO,...
                       'IN_SITU',IN_SITU);

STEP.Variables = struct('title',blockTitle,'description',description);
STEP.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
PFI.Children = struct('INFORMATION',INFORMATION,'INITIAL_GEOMETRY_DEFINITION',INITIAL_GEOMETRY_DEFINITION,'STEP',STEP);
PFI.Variables = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
nodeNames = {   'PFI'; ...
                'INFORMATION';...
                'INITIAL_GEOMETRY_DEFINITION';...
                'STEP';...
                'SOLUTION_CONTROL';...
                'LINEAR_SOLVE';...
                'TIME_CONTROL';...
                'FLOW_DEFINITION';...
                'REYNOLDS';...
                'CIRCULATION';...
                'RESULT';...
                'LOG_INFO';...
                'IN_SITU'};

nodeProperties = {  rmfield(PFI,'Children');...
                    rmfield(INFORMATION,'Children');...
                    rmfield(INITIAL_GEOMETRY_DEFINITION,'Children');...
                    rmfield(STEP,'Children');...
                    rmfield(SOLUTION_CONTROL,'Children');...
                    rmfield(LINEAR_SOLVE,'Children');...
                    rmfield(TIME_CONTROL,'Children');...
                    rmfield(FLOW_DEFINITION,'Children');...
                    rmfield(REYNOLDS,'Children');...
                    rmfield(CIRCULATION,'Children');...
                    rmfield(RESULT,'Children');...
                    rmfield(LOG_INFO,'Children');...
                    rmfield(IN_SITU,'Children')};                    
                    
                    
ParentChild = [ ...
                2 1;
                3 1;
                4 1;...
                5 4;...
                6 5;...
                7 5;...
                8 4;...
                9 8;...
                10 8;...
                11 4;...
                12 4;...
                13 4;];

nodeTable = table(nodeNames,nodeProperties,'VariableNames',{'Name','SectionProperties'});
edgeTable = table([ParentChild(:,2) ParentChild(:,1)],'VariableNames',{'EndNodes'});
pGraph = digraph(edgeTable,nodeTable);
plot(pGraph)