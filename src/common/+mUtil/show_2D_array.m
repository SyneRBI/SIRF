function show_2D_array(array, the_title, x_label, y_label)
% Displays a 2D array.

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
scale = max(abs(array(:)))/255;
image = uint8(array(:,:)/scale);
figure;
imshow(image, 'Colormap', jet(255))
xlabel(y_label)
ylabel(x_label)
set(gca, 'XTick', [1 ny])
set(gca, 'YTick', [1 nx])
axis on
title(the_title);
