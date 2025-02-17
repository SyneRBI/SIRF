try:
    import torch
except ImportError:
    raise ImportError('Failed to import torch. Please install PyTorch first.')

# based on 
# https://github.com/educating-dip/pet_deep_image_prior/blob/main/src/deep_image_prior/torch_wrapper.py


class _objectiveFunctionModule3D(torch.autograd.Function):
    @staticmethod
    def forward( ctx, x, image_template, sirf_obj):
        ctx.device = x.device
        ctx.sirf_obj = sirf_obj
        ctx.image_template = image_template
        ctx.x = x.detach().cpu().numpy().squeeze()
        ctx.x = ctx.image_template.fill(ctx.x)
        value_np = ctx.sirf_obj.get_value(ctx.x)
        return torch.tensor(value_np).to(ctx.device)

    @staticmethod
    def backward(
            ctx, 
            in_grad):
        grads_np = ctx.sirf_obj.get_gradient(ctx.x).as_array()
        grads = torch.from_numpy(grads_np).to(ctx.device) * in_grad
        return grads.unsqueeze(dim=0), None, None, None

class _AcquisitionModelForward(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, image_template, data_template, sirf_acq_model):
        ctx.sirf_acq_model = sirf_acq_model
        ctx.image_template = image_template
        ctx.data_template = data_template
        x_sirf = sirf_to_torch(x, image_template)
        x_np = x.detach().cpu().numpy()
        x_np = ctx.image_template.fill(x_np[None])
        proj_data_np = ctx.sirf_acq_model.forward(x_np).as_array()
        proj_data = torch.from_numpy(proj_data_np).requires_grad_().to(x.device)
        return proj_data

    @staticmethod
    def backward(ctx, data):
        data_np = data.detach().cpu().numpy()
        data_np = ctx.data_template.fill(data_np)
        grads_np = ctx.sirf_acq_model.backward(data_np).as_array()
        grads = torch.from_numpy(grads_np).requires_grad_().to(data.device)
        return grads, None, None, None, None

class AcquisitionFactory():
    """ I am not too sure this is a good idea...
    So the idea here is the following:
    - We have a sirf object that is the acquisition model
    - This acquisition model has the same geometry, but may have additional
    components that vary between samples from the dataset.
    We need to be able to choose the correct the acquisition model wrapper based
    on the data, also there is the assumption that we have all the same data
    components for every sample in the dataset, but heyho this is a start."""
    def __init__(self, acq_model, sirf_data_in_template, sirf_data_out_template, device, *args):
        # find a way of checking if the input or output is an image
        # then throw a warning and choose whether the wrapper should be forward
        # or backward
        if len(args) > 0:
            # detemine the modality of the data
            self.modality = args[0]
        if len(args) > 1:
            # determine the additional components of the data
            self.second_arg = args[1]
            # ... and so on
        # then begin to assign the correct acquisition model wrapper



class DatasetAcquisitionModels(torch.nn.Module):
    """Class that uses the acquisition model wrapper, to separate the geometric
    forward model from components that are data-dependent. This has the is meant
    for data of the dimensions [batch, channel, *sirf_template.shape], where
    the sirf template is either the image or the measurement template."""
    def __init__(self, acq_model, sirf_data_in_template, sirf_data_out_template, device):
        super(DatasetAcquisitionModels, self).__init__()
        # PET - multiplicative and additive factors
        # MR - coils maps
        # if acq_model is a sirf object from stir
            # if acq_model has scatter: scatter == True, else scatter == False
            # if acq_model has attn: attn == True, else attn == False
            # etc etc
            # self.acq_model = AcquisitionFactory(acq_model, sirf_image_template,
            # sirf_measurement_template, devicePET, scatter, attn, device)
        # elif acq_model is a sirf object from gadgetron
            # if acq_model has coil_sens: coil_sens == True, else coil_sens == False
            # etc etc
            # self.acq_model = AcquisitionFactory(acq_model, sirf_image_template, 
            # sirf_measurement_template, deviceMR, coil_sens, device)
        # else:
            # raise ValueError('acq_model must be a sirf object from STIR or Gadgetron')
        self.acq_model = acq_model
        self.sirf_data_in_template = sirf_data_in_template
        self.sirf_data_out_template = sirf_data_out_template
    # have a forward that can take a dynamic amount of arguments
    def forward(self, data_in, *args):
        # check if data of form [batch, channel, *sirf_data_in_template.shape]
        if data.shape[2:] != sirf_data_in_template.shape:
            raise ValueError('Data must be of the form [batch, channel, *sirf_data_in_template.shape]')

        data_out = torch.zeros(data_in.shape[0], data_in.shape[1], *sirf_data_out_template.shape)
        for i in data.shape[0]:
            for j in data.shape[1]:
                if len(args) > 0:
                    additional_data = []
                    for k in range(len(args)):
                        additional_data.append(args[k][i, j])
                data_out[i, j] = self.acq_model(data_in[i, j], *additional_data)
        return self.acq_model(*args)



# Objective function
# Acquistion model
#   With torch functionals and comparison
# Jacobian vector for STIR
# similar thing for Gadgetron
# adjointness test
# gradcheck, check traceability