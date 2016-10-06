% equivalent of demo2 using python interpreter

PYPATH = py.sys.path;
PGPATH = [getenv('SRC_PATH') '/xGadgetron/pGadgetron'];
if count(PYPATH, PGPATH) == 0
    insert(PYPATH, int32(0), PGPATH);
end

try
    % acquisitions will be read from this HDF file
    file = input('raw data file: ', 's');
    input_data = py.pGadgetron.MR_Acquisitions(file);
    
    gadgets = py.list(...
        {'RemoveROOversamplingGadget', 'SimpleReconGadgetSet'});
    recon = py.pGadgetron.ImagesReconstructor(gadgets);
    recon.set_input(input_data);
    recon.process();
    complex_images = recon.get_output();
    images = complex_images.real();
    
    n = int64(images.number());
    fprintf('Enter slice number to view it\n')
    fprintf('(a value outside the range [1 : %d] will stop this loop)\n', n)
    while (true)
        i = input('slice: ');
        if i < 1 || i > n
            break
        end
        pdata = images.image_as_array(py.int(i));
        shape = py.numpy.asarray(pdata.shape);
        data = double(py.array.array('d', py.numpy.nditer(pdata)));
        dims = int64(py.array.array('i', py.numpy.nditer(shape)));
        data = squeeze(reshape(data, dims))';
        data = data/max(max(data));
        figure(i)
        imshow(data)
    end

catch err
    % display error information
    fprintf('%s\n', err.message)
    fprintf('error id is %s\n', err.identifier)
end

