function IP =genInputParsers
    IP.INFORMATION = parseINFORMATION();
    IP.INITIAL_GEOMETRY_DEFINITION = parseINITIAL_GEOMETRY_DEFINITION();
    IP.STEP = parseSTEP();
    IP.SOLUTION_CONTROL = parseSOLUTION_CONTROL();
    IP.LINEAR_SOLVE = parseLINEAR_SOLVE();
    IP.TIME_CONTROL = parseTIME_CONTROL();
    IP.FLOW_DEFINITION = parseFLOW_DEFINITION();
    IP.REYNOLDS = parseREYNOLDS();
    IP.CIRCULATION = parseCIRCULATION();
    IP.RESULT = parseRESULT();
    IP.LOG_INFO = parseLOG_INFO();
    IP.IN_SITU = parseIN_SITU();    
end

%%
function IP = parseINFORMATION
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end

%%
function IP = parseINITIAL_GEOMETRY_DEFINITION
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addRequired(IP,'nx',@(x)validateattributes(x,'double',{'scalar','integer','>',4}))
addRequired(IP,'ny',@(x)validateattributes(x,'double',{'scalar','integer','>',4}))
addRequired(IP,'xmin',@(x)validateattributes(x,'double',{'isfinite','isreal'}))
addRequired(IP,'xmax',@(x)validateattributes(x,'double',{'isfinite','isreal'}))
addRequired(IP,'ymin',@(x)validateattributes(x,'double',{'isfinite','isreal'}))
addRequired(IP,'ymax',@(x)validateattributes(x,'double',{'isfinite','isreal'}))
end

%%
function IP = parseSTEP
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
%
addRequired(IP,'SOLUTION_CONTROL',@(x)validateattributes(x,{'struct'}))
addRequired(IP,'FLOW_DEFINITION',@(x)validateattributes(x,{'struct'}))
addOptional(IP,'RESULT',@(x)validateattributes(x,{'struct'}))
addOptional(IP,'LOG_INFO',@(x)validateattributes(x,{'struct'}))
addOptional(IP,'IN_SITU',@(x)validateattributes(x,{'struct'}))
end

%%
function IP = parseSOLUTION_CONTROL
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addParameter(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addParameter(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
%
addRequired(IP,'LINEAR_SOLVE',@(x)validateattributes(x,{'struct'}))
addRequired(IP,'TIME_CONTROL',@(x)validateattributes(x,{'struct'}))
end

%%
function IP = parseLINEAR_SOLVE
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
Solvers = {"Jacobi", "PCG","GMRES","Direct","Linfactor","Inverse"};
Hardware = {"CPU", "GPU"};
Preconditioners = {"ILU"};
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addParameter(IP,'relTol',1e-5,@(x)validateattributes(x,{'double'},{'positive'}))
addParameter(IP,'maxIter',1e3,@(x)validateattributes(x,{'double'},{'positive','integer'}))
addParameter(IP,'solver',"Jacobi",@(x)validatestring(x,Solvers))
addParameter(IP,'hardware',"CPU",@(x)validatestring(x,Hardware))

addOptional(IP,'preconditioner',@(x)validatestring(x,Preconditioners))
addOptional(IP,'kspace',@(x)validateattributes(x,{'double'},{'positive','integer'}))
end

%%
function IP = parseTIME_CONTROL
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end

%%
function IP = parseFLOW_DEFINITION
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
%
addRequired(IP,'REYNOLDS',@(x)validateattributes(x,{'struct'}))
addRequired(IP,'CIRCULATION',@(x)validateattributes(x,{'struct'}))
end

%%
function IP = parseREYNOLDS
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
funcTypes = {'constant'};
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addRequired(IP,'func',@(x)validatestring(x,{'string'},funcTypes))
addRequired(IP,'value',@(x)validateattributes(x,{'double'},{'positive','finite'}))
addParameter(IP,'scale',@(x)validateattributes(x,{'double'},{'positive','finite'}))
end

%%
function IP = parseCIRCULATION
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end

%%
function IP = parseRESULT
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end

%%
function IP = parseLOG_INFO
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end

%%
function IP = parseIN_SITU
IP = inputParser;
IP.CaseSensitive = false;
IP.KeepUnmatched = false;
IP.PartialMatching = false;
IP.StructExpand = true;
%
addOptional(IP,'title','',@(x)validateattributes(x,{'string'},{'nonempty'}))
addOptional(IP,'description','',@(x)validateattributes(x,{'string'},{'nonempty'}))
end
