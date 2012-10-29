function PFI_GUI
% The following line is primarily for ease of debugging within MATLAB
clear; close all; clc

%% Generate the figure based on the user's screen resolution/size
% Get screensize of root window... ie screen #1
scrsz=get(0,'ScreenSize'); % scrsz = [xloc yloc width height] pixels
hf=figure('MenuBar','none','NumberTitle','off','Resize','off','Position',...
    [scrsz(3)/8,scrsz(4)/6,400,600],'Color','w');
movegui(hf,'center')  % Move the figure to the center of screen

%% Generate the Text
% uicontrol text object with title's text
uicontrol('Parent',hf,'Style','text','String',[{'Welcome to'};...
    {'Parabola Flow Interactive'}],'Position',[1,520,399,80],'FontSize',...
    24,'BackgroundColor','w');

% uicontrol text object for description text
Info=uicontrol('Parent',hf,'Style','text','Position',[1,400,400,80],...
    'Fontsize',12,'BackgroundColor','w');

% Cell containing text for uicontrol text object: Info
InfoSTR={'This program solves the Vorticity-Streamline equations for'...
    'the flow of an incompressible fluid around a canonic parabola at'...
    'various modified Reynolds Numbers and Circulation paramters'};
% Generate wrapped text data
[InfoOutSTR,InfoPos]=textwrap(Info,InfoSTR);
% Set Generated wrapped text data as text for uicontrol text object: Info 
set(Info,'String',InfoOutSTR,'Position',...
    [(400-InfoPos(3))/2, 520-InfoPos(4),InfoPos(3),InfoPos(4)])

%%  Display Authors
% uicontrol text object for author names
Author=uicontrol('Parent',hf,'Style','text','Position',...
    [1,1,400,80],'Fontsize',14,'BackgroundColor','w');
% Cell containing text for uicontrol text object: Author
AuthorSTR=[{'By Wallace J. Morris II,Ph.D'};...
    {'MATLAB Version By Gregory J. Vernon'}];
% Generate wrapped text data
[AuthorOutSTR,AuthorPos]=textwrap(Author,AuthorSTR);
% Set Generated wrapped text data as text for uicontrol text object: Author
set(Author,'String',AuthorOutSTR,'Position',...
    [(400-AuthorPos(3))/2, 1,AuthorPos(3),AuthorPos(4)])

%% Simulation type 

% uicontrol pushbutton to start new simulation
uicontrol('Parent',hf,'Style','PushButton','String',...
    'Start New Simulation','FontSize',24,'Position',...
    [0,420-InfoPos(4),400,100],'CallBack',@NEWSIM);

% uicontrol pushbutton to continue from previous simulation output file
uicontrol('Parent',hf,'Style','PushButton','String',...
    'Continue Simulation','FontSize',24,'Position',...
    [0,320-InfoPos(4),400,100],'CallBack',@OLDSIM);
end

%%  Callback Functions
function NEWSIM(varargin)
    PFI(0) % Run PFI with new simulation data
end

function OLDSIM(varargin)
    PFI(1) % Run PFI with previous simulation output file
end