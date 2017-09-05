function [inGraph,fLines] = processInput(inFile)
%% Read the text input file
fLines = fileread(inFile);

% Separate the file by individual lines
fLines = string(strip(strsplit(fLines,'\n')'));

%% Convert each line to UPPERCASE
% Compatibility between 2016a and 2016b
if verLessThan('matlab','9.2')
    for ii = 1:length(fLines)
        fLines{ii} = upper(fLines{ii});
    end
else
    fLine = upper(fLines);
end

%% Create Structure
depth = 0;
sCount = 0;
for f = 1:length(fLines)
    if regexp(fLines{f},'BEGIN')
        sCount = sCount + 1;
        depth = depth + 1;
        sParts = string(strsplit(fLines{f},' '));
        sType(sCount,1).Name = sParts(end);
        sType(sCount,1).Depth = depth;
        sType(sCount,1).startLine = f;
    elseif regexp(fLines{f},'END')
        depth = depth - 1;
    end
end

%% Find each section's end lines
for ii = 1:length(sType)
    endString = ["END" + " " + sType(ii).Name];
    searchLines = fLines(sType(ii).startLine:end);
    endLines = strcmpi(endString,searchLines);
    endLine = find(endLines,1,'first');
    sType(ii).endLine = (sType(ii).startLine -1) + endLine;
end

%% Make unique node names
[uNames, uIdx] = unique([sType.Name]);
isDuplicate = ~ismember([1:length(sType)],uIdx);
dNames = unique([sType(isDuplicate).Name]);
for ii = 1:length(dNames)
    fDuplicates = find(strcmpi([sType.Name],dNames(ii)));
    for jj = 1:length(fDuplicates)
        sType(fDuplicates(jj)).Name = dNames(ii) + "-" + num2str(jj);
    end
end

%% Extract section field-value pairs
for ii = 1:length(sType)
    startLine = sType(ii).startLine;
    endLine = sType(ii).endLine;
    
    % Remove lines above section
    notSectionLines = [1:startLine];
    % Remove lines below section
    notSectionLines = [notSectionLines [endLine:length(fLines)]];
    % Remove empty lines
    notSectionLines = [notSectionLines find(fLines == "")'];
    for jj = 1:length(sType)
        if ii == jj
            continue
        else
            % Remove lines in enclosed subsections
            if sType(jj).startLine > startLine && sType(jj).endLine < endLine
                notSectionLines = [notSectionLines [sType(jj).startLine:sType(jj).endLine]];
            end
        end
    end
    SectionLines = fLines;
    SectionLines(notSectionLines) = [];
    
    if ~isempty(SectionLines)
        for p = 1:length(SectionLines)
            fvp = string(strsplit(SectionLines{p},"="));
            sType(ii).field(p) = strip(fvp(1));
            sType(ii).value(p) = strip(fvp(2));
        end
    end
end

%% Create Graph of Simulation Input Structure
% Define NodeTable
for ii = 1:length(sType)
    nodeName{ii,1} = char(sType(ii).Name);
    FVP{ii,1}.Names = sType(ii).field';
    FVP{ii,1}.Values = sType(ii).value';
end
% Define Edges
eCount = 0;
for ii = 2:length(sType)
    foundParent = false;
    jj = 0;
    while foundParent == false
        jj = jj+1;
        if sType(ii-jj).Depth < sType(ii).Depth
            foundParent = true;
            eCount = eCount+1;
            edgeNodes(eCount,:) = [ii-jj ii];
        end
    end
end
EdgeTable = table(edgeNodes,'VariableNames',{'EndNodes'});
NodeTable = table(nodeName,FVP,'VariableNames',{'Name','Parameters'});
inGraph = digraph(EdgeTable,NodeTable);
plot(inGraph)


%% Verify input parameters
inGraph = verifyInput(inGraph);

end