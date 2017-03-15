% Forward projection demo: creates an image, forward-projects it to simulate
% acquisition data and uses this data to reconstruct this image

set_up_pet
import PET.*

try
    % info() printing suppressed, warning() and error() print to stdout
    printer = Printer();
    % all printing goes to stdout 
    % printer = Printer('stdout');
    % info() prints to file
    % printer = Printer('stir_demo4_info.txt');
    % info() and warning() print to file
    % printer = Printer('stir_demo4_info.txt', 'stir_demo4_warn.txt');
    % all printing goes to files
    % printer = Printer...
    %     ('stir_demo4_info.txt', 'stir_demo4_warn.txt', 'stir_demo4_errr.txt');

    % create empty image
    image = ImageData();
    image_size = [111, 111, 31];
    voxel_size = [3, 3, 3.375];
    image.initialise(image_size, voxel_size)

    % add ellipsoidal cylinders
    shape = EllipsoidalCylinder();

    shape.set_length(400);
    shape.set_radii([40, 100]);
    shape.set_origin([60, 0, 10]);
    image.add_shape(shape, 1.0)

    shape.set_radii([30, 30])
    shape.set_origin([-30, 60, 10])
    image.add_shape(shape, 1.5)

    shape.set_origin([-30, -60, 10])
    image.add_shape(shape, 0.75)

    % z-coordinate of the xy-section to plot
    z = int32(image_size(3)/2);

    % plot the image
    data = image.as_array();
    figure(1)
    data = data/max(max(max(data)));
    imshow(data(:,:,z));
    title('phantom')

    % define the acquisition model
    am = AcquisitionModelUsingMatrix();
    % forward-project the image to obtain 'raw data';
    % raw data selected by the user is used as a template
    [filename, pathname] = uigetfile('*.hs', 'Select raw data file', pet_data_path);
    templ = AcquisitionData(fullfile(pathname, filename));
    fprintf('setting up acquisition model...\n')
    am.set_up(templ, image)
    fprintf('projecting...\n')
    fd = am.forward(image); % fd sits in memory
%     fd = am.forward(image, 'demo4data.hs'); % fd sits in this file

    data = fd.as_array();
    x = int32(image_size(1)/2);
    y = int32(image_size(2)/2);
    fprintf('data(%d,%d,%d) = %f\n', x, y, z, data(x, y, z))
    fd.fill(10*data)
    data = fd.as_array();
    fprintf('data(%d,%d,%d) = %f\n', x, y, z, data(x, y, z))
    ad = AcquisitionData(fd);
    ad.fill(fd)
    data = ad.as_array();
    figure(2)
    imshow(data(:,:,z)/max(max(max(data))));
    title('simulated acquisition data')
    
    fprintf('back-projecting...\n')
    update = am.backward(fd);
    data = update.as_array();
    figure(3)
    imshow(data(:,:,z)/max(max(max(data))));
    title('back-projection of simulated data')

catch err
    % display error information
    fprintf('??? %s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
