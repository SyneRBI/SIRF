function acquisitions_handling(engine)
% ACQUISITIONS_HANDLING Demo illustrating acquisitions pre-processing 
% and plotting.
%
% In MATLAB, there are also ISMRMRD tools available for examining 
% data before processing.
%


% Select and import SIRF MATLAB MR package so that SIRF MR objects can be 
% created in this function without using the prefix 'MR.'
%set_up_mr
%import MR.*
if nargin < 1
    engine = [];
end
eval(setup_MR(engine))

% acquisitions will be read from an HDF file
[filename, pathname] = uigetfile('*.h5', 'Select raw data file', mr_data_path);
input_data = AcquisitionData(fullfile(pathname, filename));

na = input_data.number();
fprintf('%d acquisitions (readouts) found\n', na)

fprintf('sorting acquisitions...\n')
input_data.sort()

% dimensions of k-space for one slice?
[ns, ny, nc] = input_data.slice_dimensions();

% dimensions method returns size of all (i.e. including noise data) if 
% argument is passed in or if 'all' is passed in. Passing in anything else
% means not all !!
[ns, nc, na] = input_data.dimensions('not all');

% pre-process acquisition data
fprintf('processing acquisitions...\n')
processed_data = preprocess_acquisitions(input_data);
processed_data.sort()

% selected methods for getting information
flags = input_data.get_info('flags');
es1 = input_data.get_info('encode_step_1');
slice = input_data.get_info('slice');
repetition = input_data.get_info('repetition');

while true
    num = input('enter acquisition number (0 to stop this loop): ');
    if num < 1 || num > na
        break
    end
    a = input_data.acquisition(num); % There are methods for getting selected 
                                     % properties such as flags, but not
                                     % the complete set (user may be better
                                     % using MATLAB ISMRMRD tools)
    fprintf('samples: %d\n', a.number_of_samples())
    fprintf('flags: %d %d\n', a.flags(), flags(num))
    fprintf('encode step 1: %d %d\n', a.idx_kspace_encode_step_1(), es1(num))
    fprintf('slice: %d %d\n', a.idx_slice(), slice(num))
    fprintf('repetition: %d %d\n', a.idx_repetition(), repetition(num))
end

% Data returned as complex array
iarr0 = input_data.as_array();
parr0 = processed_data.as_array();

is = ns/2;
ic = nc/2;
ia = na/2;
disp(['Value of one array element: '])
disp(iarr0(is, ic, ia))

iarr0(is, ic, ia) = iarr0(is, ic, ia)*10;

% Data can be replaced using fill method
input_data.fill(iarr0);
processed_data.fill(parr0);

iarr = input_data.as_array();
parr = processed_data.as_array();

disp(['Value of same array element after replacement with 10x data: '])
disp(iarr(is, ic, ia))

iarr = abs(iarr);
parr = abs(parr);
iarr = iarr/max(max(max(iarr)));
parr = parr/max(max(max(parr)));

nz = idivide(na, ny);

has_ipt = exist('imadjust','file') && exist('mat2gray','file');

fprintf...
    ('Enter z-coordinate of the slice to view the acquired data for it\n')
fprintf...
    ('(a value outside the range [1 : %d] will stop this loop)\n', nz)
while (true)
    z = int32(input('slice: '));
    if z < 1 || z > nz
        break
    end
    fprintf('Enter coil number to view the acquired data for it\n')
    fprintf...
        ('(a value outside the range [1 : %d] will stop this loop)\n', nc)
    while (true)
        c = input('coil: ');
        if c < 1 || c > nc
            break
        end
        figure('Name',['Input Data as array, coil: ',num2str(c)])
        data = squeeze(iarr(:, c, ny*(z - 1) + 1 : ny*z));
        if has_ipt
            imshow(imadjust(mat2gray(abs(data)),[0 0.7],[],0.2))
        else
            imshow(data)
        end
        title('Input data')
        cc = double(c + nc);
        figure('Name',['Processed data as array, coil: ',num2str(c)])
        data = squeeze(parr(:, c, ny*(z - 1) + 1 : ny*z));
        if has_ipt
            imshow(imadjust(mat2gray(abs(data)),[0 0.7],[],0.2))
        else
            imshow(data)
        end
        title('Processed data')
    end
end


