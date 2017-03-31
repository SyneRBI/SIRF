function show_3D_array(array, sup_title, label)

shape = size(array);
nx = shape(1);
ny = shape(2);
nz = int16(shape(3));
s = double(nz*nx/ny);
rows = int16(round(sqrt(s)));
if rows < 1
    rows = 1;
elseif rows > nz
    rows = nz;
end
cols = idivide(nz - 1, rows) + 1;
rows = double(rows);
cols = double(cols);
h = figure;
for z = 1 : nz
    p = double(z);
    subplot(rows, cols, p)
    imshow(array(:,:,z))
    if nargin > 2
        sub_title = sprintf('%s %d', label, z);
    else
        sub_title = sprintf('%d', z);
    end
    title(sub_title)
end
%text(-10, -10, sup_title)
h.NextPlot = 'add';
a = axes;
ht = title(sup_title);
a.Visible = 'off';
ht.Visible = 'on';
set(ht, 'FontSize', 16);
end