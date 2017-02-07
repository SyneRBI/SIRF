% Upper-level demo, GRAPPA reconstruction of undersampled data.

select_gadgetron

try

    % define raw data source
    [filename, pathname] = uigetfile('*.h5', 'Select raw data file');
    input_data = AcquisitionData(fullfile(pathname, filename));

    % pre-process acquisitions
    preprocessed_data = preprocess_acquisitions(input_data);

    % perform reconstruction
    recon = GenericCartesianGRAPPAReconstruction();
    % for undersampled acquisition data GRAPPA can compute G-factors
    % in addition to reconstructed images
%    recon.compute_gfactors(false);
    recon.compute_gfactors(true);
    recon.set_input(preprocessed_data);
    fprintf('---\n reconstructing...\n');
    recon.process();
    images = recon.get_output('image');
    gfacts = recon.get_output('gfactor');

    idata = abs(images.as_array());
    gdata = abs(gfacts.as_array());
    idata = idata/max(max(max(idata)));
    gdata = gdata/max(max(max(gdata)));

    n = images.number();
    fprintf('Enter slice number to view its data\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop).\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        figure(i)
        imshow(idata(:,:,i));
        title(['image ' num2str(i)])
        %if i <= gfacts.number()
        figure(i + n)
        imshow(gdata(:,:,i));
        title(['G-factor ' num2str(i)])
        %end
    end
    
catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end
