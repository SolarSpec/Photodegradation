clear all

%MAKE SURE FILES NAMES ARE
%Trial#"whatever you want"t=##
% example: Trial1 bulkCnx t=10

%Message box to confirm the files are all named properly
uiwait(msgbox({'Please Check files Names, MUST be of the format';'FIRST: Trial#'; 'AND CONTAIN: t=##_';'';'Example: Trial1 Bulk CNx t=20_Absorbance_15-12-04'},'modal'));
uiwait(msgbox({'If file name is not as previously specificed'; 'STOP code and fix them or this code will not work properly'},'modal'));

%List files in current working directory
files=dir('*_Absorbance_*');

%The number of files
numfile=length(files);
absorbanceData=cell(1,numfile);

%import Data
for k=1:numfile
    filename=files(k).name; %gives the name of the kth file
    absorbanceData{k}=importdata(filename); %imports that file to the kth cell of the previously empty cell array
    absorbanceData{1, k} = setfield(absorbanceData{1,k},'time',extractBetween(filename,"=","_"));
end

%Input box so you write in the number of trials in the working directory
prompt = {'Enter The Number of Trials:'};
dlgtitle = 'Input';
dims = [1 50];
numtrial = inputdlg(prompt,dlgtitle,dims);
numtrial=str2double(numtrial{1,1});
trialAbsData=cell(numtrial,numfile/numtrial) ;

%This is splitting up the files based on the name, into trial 1 trial 2
for j=1:numtrial;
    b=1;
    for L=1:numfile;
        i=files(L).name(6);
        if str2double(i)==j;
            trialAbsData{j,b}=absorbanceData{1,L};
            b=b+1;
        end   
    end
end

%find the max absorbance and coresponding wavelength, then place into a
%'final data' file where you get the data for each trial
j=0;
[c,d]= size(trialAbsData);
%For the number of Trials
for count=1:numtrial
    %for the number of data sets in each trial
    for k = 1:d
        b=0;
        % counting number of values below 400nm in each data set
       while b < 400;
           j = j + 1 ;
           b=trialAbsData{count, k}.data(j,1);
       end
       %Now find max wavelength and corresponding absorabance
       [abs,I]= max(trialAbsData{count, k}.data(j:800,2));
       %This is the row number the max absorbance is in
       value = (I + (j-1));
       Lambamaxvalues{k,(count+(count-1))} = trialAbsData{1, 1}.data(value,1);
       Lambamaxvalues{k, (count+count)} = abs;
       %Get the time of trial for the x-axis and add to final data
       timestamp =(trialAbsData{count, k}.time(1));
       finaldata{k,count+(2*(count-1))}=str2num(timestamp{1});
       %get y values for final data
       finaldata{k,((count+1)+2*(count-1))}= (Lambamaxvalues{k, count*2})/(Lambamaxvalues{1, count*2});
       finaldata{k,((count+2)+2*(count-1))}= log((Lambamaxvalues{1, count*2})/(Lambamaxvalues{k, count*2}));
    end
end
%Sort Final data so the times are in order from lowest to highest
for j=1:c 
    for k=1:d
        t=k;
        while t>1
            if finaldata{t, (1+3*(j-1))}<finaldata{t-1, (1+3*(j-1))} ;
                a = finaldata(t-1,[(1+3*(j-1)):(3*j)]);
                b = finaldata(t,[(1+3*(j-1)):(3*j)]) ;
                finaldata (t-1,[(1+3*(j-1)):(3*j)]) = b ;
                finaldata(t,[(1+3*(j-1)):(3*j)]) = a ;
                t=t-1;
            else
                t=0;
            end
        end
    end
end

%assign varibable for the x and y values
time=cell(1,c);
CC0=cell(1,c);
lnC=cell(1,c);

%covert x and y values to number varibles for cell, or string
for k=1:c
    time{k}=cell2mat(finaldata(:,k*3-2));
    CC0{k}=cell2mat(finaldata(:,k*3-1));
    lnC{k}=cell2mat(finaldata(:,k*3));
end

%Make figure 1
figure(1)
for j= 1:c
    time{1,j};
    CC0{1,j};
    plot(time{1,j},CC0{1,j});
    hold on
end
xlabel('Time (min)');
ylabel('C/C_0');

%If need adjust axis
%xlim([0 60]);
%ylim([0 100]);

%Legand lables
d = cell(1,c);
g=[];
for k=1:c
    d{1,k}=("Trial "+ k);
    g = [g , string(d{1,k})];
end

legend(g,'Location','best','fontsize',11)
legend('boxoff');

goodplot()

hold off

%Make figure 2
figure(2)
slopes = [];
for j= 1:c
    scatter(time{1,j},lnC{1,j},'x')
    hold on
    %Add linear line of best fit
    Fit = polyfit(time{1,j},lnC{1,j},1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
    f = polyval(Fit,time{1,j});
    %Getting the slope of the line of best fit
    slope = Fit(1);
    slopes= [slopes, slope];
    %plotting line of best fit
    plot(time{1,j}, f, '--r')
    hold on
end

%Legand lables
c=c*2;
e = cell(1,c);
f=[];
p=-1;

for k=1:c
    mod(k, 2);
    if mod(k,2)== 0;
        e{1,k}=("Line Best Fit Trial "+ (k-p-1));    
    else
        p=p+1;
        e{1,k}=("Trial "+ (k-p));
    end
    f = [f , string(e{1,k})];
end

legend(f,'Location','best','fontsize',11)
legend('boxoff');
goodplot()
hold off
xlabel('Time (min)');
ylabel('ln(C_0/C');

% returing slopes and STD of the slopes in command window
slopes
if numtrial >1
    S = std(slopes)
end

function goodplot()
ax = gca;
fig = gcf;
fig.Color = 'white';
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
ax.Color = 'white';
ax.LineWidth = 1.5;
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.Box = 'on';
ax.FontName = 'Arial';
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.Units = 'centimeter';

%ax.XScale = 'log' %Uncomment as needed to make scale logarithmic
%ax.YScale = 'log'

%ax.XRuler.Exponent = 0 %Change as needed for scientific notation exponent.
%ax.YRuler.Exponent = 0

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [0.03 0.01];



%% Set data properties

if isempty(findobj(fig,'Type','Line')) == 0
    plotLineData = findobj(fig,'Type','Line');
    NumberOfLines = size(plotLineData,1);
    ColorPallet = distinguishable_colors(NumberOfLines);
    for count = 1:1:NumberOfLines
        plotLineData(count).LineWidth = 2;
        plotLineData(count).Color = ColorPallet(count,:);
    end
end

if isempty(findobj(fig,'Type','Scatter')) == 0
    plotScatterData = findobj(fig,'Type','Scatter');
    NumberOfDataSets = size(plotScatterData,1);
     ColorPallet = distinguishable_colors(NumberOfDataSets);
    for count = 1:1:NumberOfDataSets
        plotScatterData(count).Marker = 'square';
        plotScatterData(count).SizeData = 40;
        plotScatterData(count).MarkerEdgeColor = ColorPallet(count,:);
        plotScatterData(count).MarkerFaceColor = plotScatterData(count).MarkerEdgeColor;
    end
end



%% Get rid of extra white space around figure

fig.WindowStyle = 'normal';
ax.OuterPosition = [0 0 8.5 7.5]; %Change as needed.
CurrentFigPos = fig.Position;
fig.Position = [CurrentFigPos(1) CurrentFigPos(2) 8.5 7.5]; %Change width and height as needed. 8.5 x 7.5 cm is suitable for 1 plot.

end
%% Function script to choose colors

function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end