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
%         idata = abs(input_data.slice_as_array(z));
%         pdata = abs(processed_data.slice_as_array(z));
        fprintf('Enter coil number to view the acquired data for it\n')
        fprintf...
            ('(a value outside the range [1 : %d] will stop this loop)\n', nc)
        while (true)
            c = input('coil: ');
            if c < 1 || c > nc
                break
            end
%             idata = idata/max(max(max(idata)));
%             pdata = pdata/max(max(max(pdata)));
            figure(c)
            data = squeeze(iarr(:, c, ny*(z - 1) + 1 : ny*z));
            imshow(data)
            %imshow(idata(:,:,c));
            title('oversampled data')
            cc = double(c + nc);
            figure(cc)
            data = squeeze(parr(:, c, ny*(z - 1) + 1 : ny*z));
            imshow(data)
            %imshow(pdata(:,:,c));
            title('de-oversampled data')
        end
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
