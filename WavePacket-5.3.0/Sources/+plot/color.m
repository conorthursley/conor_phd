%------------------------------------------------------------------------------
%
% View two functions in a single curve:
% Draw a colored curve of width w, specified by x and y, with 
% color mapping c which must be inside interval [0,1].
% Base color c0 (RGB) is used where abs(y) is below threshold.
% Plot curve with ordinate values horizontally offset.
% If swap is set to 1, abscissa and ordinate are swapped.
%
% HACK: We overloaded swap. If set to 2, we have a polar plot.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2009 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% see the README file for license details.

function color(x,y,c,c0,w,offset,swap)

% Check input parameters for consistency
if ~( isnumeric(x) & isnumeric(y) & isnumeric(c) )
    util.error('Part of the input is not numeric!')    
end

if ~isreal(x)
    util.error('Abscissa values have to be real!')
end

if ~isreal(y)
    util.error('Ordinate values have to be real!')
end

if ~isreal(c)
    util.error('Color values have to be real!')
end

if length(x)==length(y)
    n = length(x);
else
    util.error('Dimension mismatch!')
end

if  ~(length(y)==length(c))
    util.error('Dimension mismatch!')
end

% Using hue-saturation-value (HSV) scheme
[hue0 sat0 val0] = rgb2hsv (c0);
hue = hue0 * ones (n,1);  % arbitrary hue
sat = sat0 * ones (n,1);  % fully saturated
val = val0 * ones (n,1);  % fully dark

% Find indices where density is below threshold
if swap ~= 2
    epsilon = max(y)/1000;                                        
    not_too_small = find(y>epsilon)';
else
    not_too_small = 1:length(y);
end

% Hue must be inside interval [0,1] 
hue (not_too_small) = c(not_too_small) ;
val (not_too_small) = 1;

% Convert to red-green-blue (RGB) scheme 
rgb = hsv2rgb ([hue,sat,val]);

% Horizontal offset and swap abscissa and ordinate
if swap == 1
    dummy = y;
    y = x;
    x = dummy + offset;
else
    y = y + offset;
end

% Plot the first data point and append colored line-pieces one by one
plot(x(1),y(1))
for k=1:(n-1)
    line([x(k),x(k+1)],[y(k),y(k+1)],'LineWidth',w,'Color',rgb(k,:))
end  

