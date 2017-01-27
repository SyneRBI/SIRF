% Upper-level interface demo, illustrates acquisitions pre-processing and plotting.

set_up_mr

try
    % acquisitions will be read from this HDF file
    input_data = AcquisitionData('simulated_MR_2D_cartesian.h5');
    
    na = input_data.number();
    fprintf('%d acquisitions found\n', na)

    fprintf('sorting acquisitions...\n')
    input_data.sort()

    [ns, ny, nc] = input_data.slice_dimensions();

    [ns, nc, na] = input_data.dimensions('image');
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    processed_data = preprocess_acquisitions(input_data);
    processed_data.sort()
    
    flags = input_data.get_info('flags');
    es1 = input_data.get_info('encode_step_1');
    slice = input_data.get_info('slice');
    repetition = input_data.get_info('repetition');
    
    while true
        num = input('enter acquisition number: ');
        if num < 1 || num > na
            break
        end
        a = input_data.acquisition(num);
        fprintf('samples: %d\n', a.number_of_samples())
        fprintf('flags: %d %d\n', a.flags(), flags(num))
        fprintf('encode step 1: %d %d\n', a.idx_kspace_encode_step_1(), es1(num))
        fprintf('slice: %d %d\n', a.idx_slice(), slice(num))
        fprintf('repetition: %d %d\n', a.idx_repetition(), repetition(num))
    end
    
    iarr0 = input_data.as_array();
    parr0 = processed_data.as_array();
    
    is = ns/2;
    ic = nc/2;
    ia = na/2;
    disp(iarr0(is, ic, ia))
    iarr0(is, ic, ia) = iarr0(is, ic, ia)*10;
    
    input_data.fill(iarr0);
    processed_data.fill(parr0);

    iarr = input_data.as_array();
    parr = processed_data.as_array();

    disp(iarr(is, ic, ia))

    iarr = abs(iarr);
    parr = abs(parr);
    iarr = iarr/max(max(max(iarr)));
    parr = parr/max(max(max(parr)));

    nz = idivide(na, ny);
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
            figure(c)
            data = squeeze(iarr(:, c, ny*(z - 1) + 1 : ny*z));
            imshow(data)
            title('oversampled data')
            cc = double(c + nc);
            figure(cc)
            data = squeeze(parr(:, c, ny*(z - 1) + 1 : ny*z));
            imshow(data)
            title('de-oversampled data')
        end
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
