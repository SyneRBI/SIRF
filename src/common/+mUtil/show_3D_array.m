function show_3D_array(array, sup_title, label)

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