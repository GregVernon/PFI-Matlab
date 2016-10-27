function input = processInput(inFile)

%% Read the text input file
fLines = fileread(inFile);

% Separate the file by individual lines
fLines = strsplit(fLines,'\n');

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
    sectionType{1,n1} = parseLine{end};
  elseif regexp(fLines{ii},'END')
    n2 = n2+1;
    parseLine = strsplit(fLines{ii});
    sectionType{1,n2} = parseLine{end};
  end
end
