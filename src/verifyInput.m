function inGraph = verifyInput(inGraph)
%% Get list of legal parameters
pfiGraph = generateLegalParameters;

%% Verify PFI Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'PFI');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'PFI'));

if (isUsed == false)
    error('"PFI" block is required');
end

%% Verify INFORMATION Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'INFORMATION');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'INFORMATION'));

if (isRequired == true) && (isUsed == false)
    error('"INFORMATION" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'INFORMATION');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    end
end

%% Verify INITIAL_GEOMETRY_DEFINITION Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'INITIAL_GEOMETRY_DEFINITION');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'INITIAL_GEOMETRY_DEFINITION'));

if (isRequired == true) && (isUsed == false)
    error('"INITIAL_GEOMETRY_DEFINITION" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'INITIAL_GEOMETRY_DEFINITION');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'nx')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'ny')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'xmin')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'xmax')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'ymin')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'ymax')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    end
end

%% Verify STEP Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'STEP');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'STEP'));

if (isRequired == true) && (isUsed == false)
    error('"STEP" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'STEP');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    end
end

%% Verify SOLUTION_CONTROL Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'SOLUTION_CONTROL');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'SOLUTION_CONTROL'));

if (isRequired == true) && (isUsed == false)
    error('"SOLUTION_CONTROL" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'SOLUTION_CONTROL');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    end
end

%% Verify LINEAR_SOLVE Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'LINEAR_SOLVE');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'LINEAR_SOLVE'));

if (isRequired == true) && (isUsed == false)
    error('"LINEAR_SOLVE" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'LINEAR_SOLVE');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'relative tolerance')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'solver')
        assert(any(strcmpi(pParam.solver.validValues,iParam.Values{ii})))
    elseif strcmpi(iParam.Names(ii),'preconditioner')
        assert(any(strcmpi(pParam.precond.validValues,iParam.Values{ii})))
    elseif strcmpi(iParam.Names(ii),'k-space')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'maximum iterations')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    end
end

%% Verify TIME_CONTROL Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'TIME_CONTROL');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'TIME_CONTROL'));

if (isRequired == true) && (isUsed == false)
    error('"TIME_CONTROL" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'TIME_CONTROL');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'step time')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'max time control')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'max time scale factor')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    end
end

%% Verify FLOW_DEFINITION Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'FLOW_DEFINITION');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'FLOW_DEFINITION'));

if (isRequired == true) && (isUsed == false)
    error('"FLOW_DEFINITION" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'FLOW_DEFINITION');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    end
end

%% Verify REYNOLDS Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'REYNOLDS');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'REYNOLDS'));

if (isRequired == true) && (isUsed == false)
    error('"REYNOLDS" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'REYNOLDS');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'function')
        assert(any(strcmpi(pParam.precond.validValues,iParam.Values{ii})))
    elseif strcmpi(iParam.Names(ii),'value')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'scale')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    end
end

%% Verify CIRCULATION Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'CIRCULATION');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'CIRCULATION'));

if (isRequired == true) && (isUsed == false)
    error('"CIRCULATION" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'CIRCULATION');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'function')
        assert(any(strcmpi(pParam.precond.validValues,iParam.Values{ii})))
    elseif strcmpi(iParam.Names(ii),'value')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'scale')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    end
end

%% Verify RESULT Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'RESULT');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(contains(inGraph.Nodes{:,'Name'},'RESULT'));

if (isRequired == true) && (isUsed == false)
    error('"RESULT" block is required');
else
    iRow = find(contains(inGraph.Nodes{:,'Name'},'RESULT'));
    iParam = inGraph.Nodes{iRow,'Parameters'};
end

pFieldnames = string(fieldnames(pParam));
for idx = 1:length(iRow)
    for ii = 1:length(iParam{idx}.Names)
        ipMatch = false(size(pFieldnames));
        for jj = 1:length(pFieldnames)
            ipMatch(jj) = strcmpi(iParam{idx}.Names(ii), pParam.(pFieldnames{jj}).name);
        end
        
        if any(ipMatch) == false
            error(['"' char(iParam{idx}.Names(ii)) '" not a valid command'])
        end
        
        inGraph.Nodes{iRow(idx),'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
        if strcmpi(iParam{idx}.Names(ii),'title')
            assert(isstring(iParam{idx}.Values{ii}))
        elseif strcmpi(iParam{idx}.Names(ii),'description')
            assert(isstring(iParam{idx}.Values{ii}))
        elseif strcmpi(iParam{idx}.Names(ii),'output filename')
            assert(isstring(iParam{idx}.Values{ii}))
        elseif strcmpi(iParam{idx}.Names(ii),'active step')
            assert(isstring(iParam{idx}.Values{ii}))
        elseif strcmpi(iParam{idx}.Names(ii),'output initial step')
            assert(any(strcmpi(["true","false"],iParam{idx}.Values{ii})))
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2num(char(iParam{idx}.Values{ii}));
        elseif strcmpi(iParam{idx}.Names(ii),'output final step')
            assert(any(strcmpi(["true","false"],iParam{idx}.Values{ii})))
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2num(char(iParam{idx}.Values{ii}));
        elseif strcmpi(iParam{idx}.Names(ii),'time increment')
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2double(iParam{idx}.Values{ii});
        elseif strcmpi(iParam{idx}.Names(ii),'buffer increments')
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2double(iParam{idx}.Values{ii});
        elseif strcmpi(iParam{idx}.Names(ii),'variable')
            assert(isstring(iParam{idx}.Values{ii}))
        elseif strcmpi(iParam{idx}.Names(ii),'iX')
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2num(char(iParam{idx}.Values{ii}));
        elseif strcmpi(iParam{idx}.Names(ii),'iY')
            inGraph.Nodes{iRow(idx),'Parameters'}{1}.Values{ii} = str2num(char(iParam{idx}.Values{ii}));
        end
    end
end

%% Verify LOG_INFO Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'LOG_INFO');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'LOG_INFO'));

if (isRequired == true) && (isUsed == false)
    error('"LOG_INFO" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'LOG_INFO');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'log filename')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'time increment')
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2double(iParam.Values{ii});
    elseif strcmpi(iParam.Names(ii),'initial step')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
    elseif strcmpi(iParam.Names(ii),'final step')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
    elseif strcmpi(iParam.Names(ii),'log data')
        assert(isstring(iParam.Values{ii}))
    end
end

%% Verify IN_SITU Block
pRow = strcmpi(pfiGraph.Nodes{:,'Name'},'IN_SITU');
pParam = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.Variables;
isRequired = pfiGraph.Nodes{pRow,'SectionProperties'}{1}.isRequired;
isUsed = any(strcmpi(inGraph.Nodes{:,'Name'},'IN_SITU'));

if (isRequired == true) && (isUsed == false)
    error('"IN_SITU" block is required');
else
    iRow = strcmpi(inGraph.Nodes{:,'Name'},'IN_SITU');
    iParam = inGraph.Nodes{iRow,'Parameters'}{1};
end

pFieldnames = string(fieldnames(pParam));
for ii = 1:length(iParam.Names)
    ipMatch = false(size(pFieldnames));
    for jj = 1:length(pFieldnames)
        ipMatch(jj) = strcmpi(iParam.Names(ii), pParam.(pFieldnames{jj}).name);
    end
    
    if any(ipMatch) == false
        error(['"' char(iParam.Names(ii)) '" not a valid command'])
    end
    
    inGraph.Nodes{iRow,'Parameters'}{1}.Fields(ii,1) = string(pFieldnames{ipMatch});
    if strcmpi(iParam.Names(ii),'title')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'description')
        assert(isstring(iParam.Values{ii}))
    elseif strcmpi(iParam.Names(ii),'plotFlow')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
    elseif strcmpi(iParam.Names(ii),'plotState')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
    elseif strcmpi(iParam.Names(ii),'saveFlow')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
    elseif strcmpi(iParam.Names(ii),'saveState')
        assert(any(strcmpi(["true","false"],iParam.Values{ii})))
        inGraph.Nodes{iRow,'Parameters'}{1}.Values{ii} = str2num(char(iParam.Values{ii}));
        
    end
end
