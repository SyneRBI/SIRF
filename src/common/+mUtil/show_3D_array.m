function show_3D_array(array, sup_title, label)
% Shows 3D array as a set of xy-slices, optionally labelled.

% CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
% Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
% Copyright 2015 - 2017 University College London.
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
rows = uint16(round(sqrt(s)));
% extra_row = 0;
% if nz == 1
%     rows = uint16(2);
%     extra_row = 1;
% end
if rows < 1
    rows = uint16(1);
end
if rows > nz
    rows = uint16(nz);
end
cols = idivide(nz - 1, rows) + 1;
rows = double(rows);
cols = double(cols);
h = figure;
scale = max(abs(array(:)))/255;
for z = 1 : nz
    p = double(z); % + extra_row*cols;
    subplot(rows, cols, p)
    image = uint8(array(:,:,z)/scale);
    imshow(image, 'Colormap', jet(255))
    if nz > 1
        if nargin > 2
            sub_title = sprintf('%s %d', label, z);
        else
            sub_title = sprintf('%d', z);
        end
        %t = title(sub_title);
        t = title('');
        set(t, 'FontSize', 8);
        tpos = round(get(t, 'Position'));
        txt = text(tpos(1), tpos(2) + nx + 50, sub_title, ...
            'HorizontalAlignment', 'center');
    end
end
h.NextPlot = 'add';
a = axes;
ht = title(sup_title);
a.Visible = 'off';
ht.Visible = 'on';
%set(ht, 'FontSize', 12);
end