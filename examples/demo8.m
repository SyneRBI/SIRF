% Upper-level interface demo, illustrates acquisitions pre-processing and plotting.

if ~libisloaded('mgadgetron')
    loadlibrary('mgadgetron')
end
if ~libisloaded('mutilities')
    loadlibrary('mutilities')
end

try
    % acquisitions will be read from this HDF file
    input_data = gadgetron.MR_Acquisitions('testdata.h5');
    
    na = input_data.number();
    fprintf('%d acquisitions found\n', na)

    fprintf('sorting acquisitions...\n')
    input_data.sort()

    [ns, ny, nc] = input_data.slice_dimensions();
    
    % pre-process acquisition data
    fprintf('processing acquisitions...\n')
    processed_data = gadgetron.MR_remove_x_oversampling(input_data);

    nz = idivide(na,ny);
    fprintf...
        ('Enter z-coordinate of the slice to view the acquired data for it\n')
    fprintf...
        ('(a value outside the range [1 : %d] will stop this loop)\n', nz)
    while (true)
        z = int32(input('slice: '));
        if z < 1 | z > nz
            break
        end
        idata = abs(input_data.slice_as_array(z));
        pdata = abs(processed_data.slice_as_array(z));
        fprintf('Enter coil number to view the acquired data for it\n')
        fprintf...
            ('(a value outside the range [1 : %d] will stop this loop)\n', nc)
        while (true)
            c = input('coil: ');
            if c < 1 | c > nc
                break
            end
            idata = idata/max(max(max(idata)));
            pdata = pdata/max(max(max(pdata)));
            figure(c)
            imshow(idata(:,:,c));
            title('oversampled data')
            cc = double(c + nc);
            figure(cc)
            imshow(pdata(:,:,c));
            title('de-oversampled data')
        end
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
