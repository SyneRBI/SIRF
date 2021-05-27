function err = show_3D_array(array, sup_title, x_label, y_label, label, index)
% Shows 3D array as a set of xy-slices, optionally labelled.
%    array     : 3D array
%    suptitle  : figure title
%    xlabel    : label for x axis
%    ylabel    : label for y axis
%    label     : tile title prefix
%    index     : z-slices index, either integer array or string of the form
%              : 'a, b-c, ...', where 'b-c' is decoded as 'b, b+1, ..., c';
%              : out-of-range index value causes error (non-zero) return

% SyneRBI Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
% (http://www.ccpsynerbi.ac.uk/).
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

err = 0;

shape = size(array);
nx = shape(1);
ny = shape(2);

if ndims(array) < 3
    nz = uint16(1);
else
    nz = uint16(shape(3));
end
if nargin < 6
    n = nz;
    index = 1 : nz;
else
    if ischar(index)
        try
            index = sirf.Utilities.str_to_int_list(index);
        catch
            return
        end
    end
    n = length(index);
    if n < 1
        return
    end
    for i = 1 : n
        j = index(i);
        if j < 1 || j > nz
            err = i;
            %fprintf('z-index %d is out of range, aborting display\n', j)
            return
        end
    end
end
nz = n;
s = double(nz*ny/nx);
irows = uint16(round(sqrt(s)));
if irows < 1
    irows = uint16(1);
end
if irows > nz
    irows = uint16(nz);
end
icols = idivide(nz - 1, irows) + 1;
irows = idivide(nz - 1, icols) + 1;
rows = double(irows);
cols = double(icols);

f = figure;
panel = uipanel('Parent', f, 'BorderType', 'none');
panel.Title = sup_title;
panel.TitlePosition = 'centertop'; 
panel.FontSize = 12;
panel.FontWeight = 'bold';

vmin = min(array(:));
vmax = max(array(:));
for z = 1 : nz
    row = idivide(z - 1, icols) + 1;
    col = z - (row - 1)*icols;
    p = double(z);
    subplot(rows, cols, p, 'Parent', panel)
    imshow(array(:,:,index(z)), [vmin vmax], 'Colormap', jet(255))
    if row == irows && col == 1
        xlabel(y_label)
        ylabel(x_label)
        set(gca, 'XTick', [1 ny])
        set(gca, 'YTick', [1 nx])
        axis on
    end
    if nz > 1
        if nargin > 4
            sub_title = sprintf('%s %d', label, index(z));
            title(sub_title);
        end
    end
end
%colorbar
end