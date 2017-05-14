function [input,fLines] = processInput(inFile)
%% Get list of legal parameters
PARAM = generateLegalParameters;
%% Read the text input file
fLines = fileread(inFile);

% Separate the file by individual lines
fLines = deblank(strsplit(fLines,'\n')');

%% Convert each line to UPPERCASE
% Compatibility between 2016a and 2016b
if verLessThan('matlab','9.2')
  for ii = 1:length(fLines)
    fLines{ii} = upper(fLines{ii});
  end
else
  fLine = upper(fLines);
end

%% Find "BEGIN" and "END" sections
n1 = 0;
n2 = 0;
for ii = 1:length(fLines)
  if regexp(fLines{ii},'BEGIN')
    n1 = n1+1;
    parseLine = strsplit(fLines{ii});
    sectionType{1,n1,1} = parseLine{end};
    sectionType{1,n1,2} = ii;
  elseif regexp(fLines{ii},'END')
    n2 = n2+1;
    parseLine = strsplit(fLines{ii});
    sectionType{2,n2,1} = parseLine{end};
    sectionType{2,n2,2} = ii;
  end
end

input = [];
