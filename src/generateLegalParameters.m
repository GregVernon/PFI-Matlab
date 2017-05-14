function PARAM = generateLegalParameters
%%
simTitle = struct('name','title','class','char','isRequired',false);

INFORMATION.Children = [];
INFORMATION.Variables = struct('simTitle',simTitle);
INFORMATION.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
nx = struct('name','nx','class',{'double'},'isRequired',true);
ny = struct('name','ny','class',{'double'},'isRequired',true);

xmin = struct('name','xmin','class',{'double'},'isRequired',true);
xmax = struct('name','xmax','class',{'double'},'isRequired',true);
ymin = struct('name','ymin','class',{'double'},'isRequired',true);
ymax = struct('name','ymax','class',{'double'},'isRequired',true);
GEOMETRY_DEFINITION.Children = [];
GEOMETRY_DEFINITION.Variables = struct('nx',nx,'ny',ny,'xmin',xmin,'xmax',xmax,'ymin',ymin,'ymax',ymax);
GEOMETRY_DEFINITION.isRequired = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
relTol = struct('name','relative tolerance','class',{'double'},'isRequired',false,'default',1e-5);
solver = struct('name','solver','class','char','isRequired',false,'default','Jacobi','options',{'Jacobi';'PCG';'GMRES';'Direct';'Linfactor';'Inverse'});
precond = struct('name','preconditioner','class','char','isRequired',false,'default',[]);
kspace = struct('name','k-space','class','double','isRequired',false,'default',20);
maxit = struct('name','maximum iterations','class','double','isRequired','false','default',1e4);

LINEAR_SOLVE.Children = [];
LINEAR_SOLVE.Variables = struct('relTol',relTol,'solver',solver,'precond',precond,'kspace',kspace,'maxit',maxit);
LINEAR_SOLVE.isRequired = true;
%%
stepTime = struct('name','step time','class','double','isRequired',true);
maxTime = struct('name','max time control','class','char','isRequired',false,'default','Reynolds');
scale = struct('name','max time scale factor','class','double','isRequired',false,'default',0.01);

TIME_CONTROL.Children = [];
TIME_CONTROL.Variables = struct('stepTime',stepTime,'maxTime',maxTime,'scale',scale);
TIME_CONTROL.isRequired = true;
%%
SOLUTION_CONTROL.Children = struct('LINEAR_SOLVE',LINEAR_SOLVE,'TIME_CONTROL',TIME_CONTROL);
SOLUTION_CONTROL.Variables = [];
SOLUTION_CONTROL.isRequired = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
func = struct('name','function','class','char','isRequired',true,'options',{'constant';},'default',{'constant'});
value = struct('name','value','class','double','isRequired',true);
scale = struct('name','scale factor','class','double','isRequired',false,'default',1);

REYNOLDS.Children = [];
REYNOLDS.Variables = struct('func',func,'value',value,'scale',scale);
REYNOLDS.isRequired = true;
%%
func = struct('name','function','class','char','isRequired',true,'options',{'constant';},'default',{'constant'});
value = struct('name','value','class','double','isRequired',true);
scale = struct('name','scale factor','class','double','isRequired',false,'default',1);

CIRCULATION.Children = [];
CIRCULATION.Variables = struct('func',func,'value',value,'scale',scale);
CIRCULATION.isRequired = true;
%%
FLOW_DEFINITION.Children = struct('REYNOLDS',REYNOLDS,'CIRCULATION',CIRCULATION);
FLOW_DEFINITION.Variables = [];
FLOW_DEFINITION.isRequired = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
outputName = struct('name','output filename','class','char','isRequired',true);
description = struct('name','description','class','char','isRequired',false);
activeStep = struct('name','active step','class',{'double';'char'},'isRequired',false,'default','all');
outInit = struct('name','output initial step','class','logical','isRequired',false,'default',true);
outFin = struct('name','output final step','class','logical','isRequired',false,'default',true);
tIncr = struct('name','time increment','class','double','isRequired',true);
dBuffer = struct('name','buffer increments','class','double','isRequired',false,'default',100);
variable = struct('name','variable','class','char','isRequired',true,'options',{'OMEGA','PSI','pU','pV'});
iX = struct('name','iX','class','double','isRequired',true);
iY = struct('name','iY','class','double','isRequired',true);

RESULT.Children = [];
RESULT.Variables = struct('outputName',outputName,'description',description,...
                          'activeStep',activeStep,'outInit',outInit,...
                          'outFin',outFin,'tIncr',tIncr,'dBuffer',dBuffer,...
                          'variable',variable,'iX',iX,'iY',iY);
RESULT.isRequired = false;

%%
logName = struct('name','log filename','class','char','isRequired',false,'default',strcat('PFI_log_',datestr(now,'mm_dd_yyyy__HH:MM:SS')));
tIncr = struct('name','time increment','class','double','isRequired',false,'default',1e3);
initStep = struct('name','initial step','class','logical','isRequired',false,'default',true);
finStep = struct('name','final step','class','logical','isRequired',false,'default',true);
logData = struct('name','log data','class','cell','isRequired',false,'default',{'Timestep';'Time';'Time Elapsed'});

LOG_INFO.Children = [];
LOG_INFO.Variables = struct('logName',logName,'tIncr',tIncr,...
                            'initStep',initStep,'finStep',finStep,...
                            'logData',logData);
LOG_INFO.isRequired = false;

%%
plotFlow = struct('name','plotFlow','class','logical','isRequired',false,'default',false);
plotState = struct('name','plotState','class','logical','isRequired',false,'default',false);
saveFlow = struct('name','saveFlow','class','logical','isRequired',false,'default',false);
saveState = struct('name','saveState','class','logical','isRequired',false,'default',false);

IN_SITU.Children = [];
IN_SITU.Variables = struct('plotFlow',plotFlow,'plotState',plotState,'saveFlow',saveFlow,'saveState',saveState);
IN_SITR.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
stepName = struct('name','title','class','char','isRequired',false);

STEP.Children = struct('SOLUTION_CONTROL',SOLUTION_CONTROL,...
                       'FLOW_DEFINITION',FLOW_DEFINITION,...
                       'RESULT',RESULT,'LOG_INFO',LOG_INFO,...
                       'IN_SITU',IN_SITU);

STEP.Variables = struct('stepName',stepName);
STEP.isRequired = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
PARAM.Children = struct('INFORMATION',INFORMATION,'GEOMETRY_DEFINITION',GEOMETRY_DEFINITION,'STEP',STEP);
PARAM.Variables = [];