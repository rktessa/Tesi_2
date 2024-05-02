clc;
clear all; 
close all;
%%
vettore = linspace(1,10,10);
matrice = zeros(10,10,10,10,10);
matrice(1,1,1,1,:) = transpose(vettore)




%%
% zeta = .5;                           % Damping Ratio
% wn = 2;                              % Natural Frequency
% sys = tf(wn^2,[1,2*zeta*wn,wn^2]); 
% 
% 
% f = figure;
% ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
% h = stepplot(ax,sys);
% setoptions(h,'XLim',[0,10],'YLim',[0,2]);
% 
% %Slider and its option 
% b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
%               'value',zeta, 'min',0, 'max',1);
% bgcolor = f.Color;
% bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
%                 'String','0','BackgroundColor',bgcolor);
% bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
%                 'String','1','BackgroundColor',bgcolor);
% bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
%                 'String','Damping Ratio','BackgroundColor',bgcolor);
% 
% 
% 
% b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2]));


%

ii = 1;                           % Damping Ratio
wn = 2;                              % Natural Frequency
x_i = linspace(1,10,10);
y_i1 = linspace(1,10,10);
y_i2 = linspace(1,20,10);
y_i3 = linspace(1,30,10);
y_i4 = linspace(1,40,10);
y_i = [y_i1;y_i2;y_i3;y_i4];




f = figure;
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
h = plot(ax,y_i(1,:));
xlim([1 10])
ylim([1 40])
% setoptions(h,'XLim',[0,10],'YLim',[0,40]);


%Slider and its option 
b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',ii, 'min',1, 'max',4, 'Value',1, ...
              'Units', 'Normalized',...
              'Callback', @print_val,...
              'SliderStep', [1/3 1]); %to make discretes steps
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','4','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','numero i','BackgroundColor',bgcolor);


c = uicontrol('Parent',f,'Style','slider','Position',[81,84,419,23],...
              'value',ii, 'min',1, 'max',4, 'Value',1, ...
              'Units', 'Normalized',...
              'Callback', @print_val,...
              'SliderStep', [1/3 1]); %to make discretes steps
cgcolor = f.Color;
cl1 = uicontrol('Parent',f,'Style','text','Position',[50,84,23,23],...
                'String','1','BackgroundColor',bgcolor);
cl2 = uicontrol('Parent',f,'Style','text','Position',[500,84,23,23],...
                'String','4','BackgroundColor',bgcolor);
cl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','numero i2','BackgroundColor',bgcolor);


% b.Callback = @(es,ed) updateSystem(h,tf(wn^2,[1,2*(es.Value)*wn,wn^2]));
% b.Callback = @(es,ed) plot(ax,y_i(es.Value,:));

function print_val(hObject,callbackdata)
    newval = hObject.Value;                         %get value from the slider
    newval = round(newval);                         %round off this value
    set(hObject, 'Value', newval);                  %set slider position to rounded off value
    disp(['Slider moved to ' num2str(newval)]);     %display the value pointed by slider
    
    % f = figure;
    % ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
    x_i = linspace(1,10,10);
    y_i1 = linspace(1,10,10);
    y_i2 = linspace(1,20,10);
    y_i3 = linspace(1,30,10);
    y_i4 = linspace(1,40,10);
    y_i = [y_i1;y_i2;y_i3;y_i4];

    plot(x_i,y_i(newval,:));
    xlim([1 10])
    ylim([1 40])
        
    % plot(ax,y_i(newval,:)
end



