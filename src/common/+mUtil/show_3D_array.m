function show_3D_array(array, sup_title, x_label, y_label, label)
% Shows 3D array as a set of xy-slices, optionally labelled.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% 
% This is software developed for the Collaborative Computational
% Project in Positron Emission Tomography and Magnetic Resonance imaging
% (http://www.ccppetmr.ac.uk/).
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

shape = size(array);
nx = shape(1);
ny = shape(2);

if ndims(array) < 3
    nz = uint16(1);
else
    nz = uint16(shape(3));
end
s = double(nz*ny/nx);
irows = uint16(round(sqrt(s)));
if irows < 1
    irows = uint16(1);
end
if irows > nz
    irows = uint16(nz);
end
icols = idivide(nz - 1, irows) + 1;
rows = double(irows);
cols = double(icols);

f = figure;
panel = uipanel('Parent', f, 'BorderType', 'none');
panel.Title = sup_title;
panel.TitlePosition = 'centertop'; 
panel.FontSize = 12;
panel.FontWeight = 'bold';

scale = max(abs(array(:)))/255;
for z = 1 : nz
    row = idivide(z - 1, icols) + 1;
    col = z - (row - 1)*icols;
    p = double(z);
    subplot(rows, cols, p, 'Parent', panel)
    if z <= nz
        image = uint8(array(:,:,z)/scale);
        imshow(image, 'Colormap', jet(255))
    else
        imshow(ones(nx, ny), 'Colormap', jet(255))
    end
    if row == irows && col == 1
        xlabel(y_label)
        ylabel(x_label)
        set(gca, 'XTick', [1 ny])
        set(gca, 'YTick', [1 nx])
        axis on
    end
    if nz > 1
        if nargin > 4
            sub_title = sprintf('%s %d', label, z);
        else
            sub_title = sprintf('%d', z);
        end
        title(sub_title);
    end
end
end