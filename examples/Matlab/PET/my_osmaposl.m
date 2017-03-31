function image = my_osmaposl...
    (init_image, obj_fun, prior, filter, num_subsets, num_subiterations)

image = init_image;
for iter = 1 : num_subiterations
    fprintf('\n------------- Subiteration %d\n', iter) 

    % select subset
    subset = iter - 1;

    % get sensitivity as ImageData
    sens_image = obj_fun.get_subset_sensitivity(subset);

    % get backprojection of the ratio of measured to estimated
    % acquisition data
    grad_image = obj_fun.get_backprojection_of_acquisition_ratio...
                 (image, subset);

    % get gradient of prior as ImageData
    prior_grad_image = prior.get_gradient(image);

    % copy to Python arrays
    image_array = image.as_array();
    sens_array = sens_image.as_array();
    grad_array = grad_image.as_array();
    prior_grad_array = prior_grad_image.as_array();

    % update image data
    denom = sens_array + prior_grad_array./num_subsets;
    eps = 1e-6*max(abs(denom(:)));
    denom(denom < eps) = eps; % avoid division by zero
    update = grad_array./denom;
    image_array = image_array.*update;

    % fill current image with new values
    image.fill(image_array);

    % apply filter
    filter.apply(image);

end
end