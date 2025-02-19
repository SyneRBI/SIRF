try:
    import torch
except ModuleNotFoundError:
    raise ModuleNotFoundError('Failed to import torch. Please install PyTorch first.')


# based on 
# https://github.com/educating-dip/pet_deep_image_prior/blob/main/src/deep_image_prior/torch_wrapper.py


def sirf_to_torch(sirf_src, device, requires_grad=False):
    if requires_grad:
        # use torch.tensor to infer data type
        return torch.tensor(sirf_src.as_array(), requires_grad=True).to(device)
    else:
        return torch.tensor(sirf_src.as_array()).to(device)

def torch_to_sirf(torch_src, sirf_dest): 
    return sirf_dest.fill(torch_src.detach().cpu().numpy())


class _ObjectiveFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            sirf_image_template, 
            sirf_obj_func
            ):

        device = torch_image.device
        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        value_np = sirf_obj_func.get_value(sirf_image_template).as_array()
        if torch_image.requires_grad:
            ctx.save_for_backward(device, sirf_image_template, sirf_obj_func)
            return torch.tensor(value_np, requires_grad=True).to(device)
        else:
            return torch.tensor(value_np).to(torch_image.device)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        device, sirf_image_template, sirf_obj_func = ctx.saved_tensors
        tmp_grad = sirf_obj.get_gradient(sirf_image_template)
        grad = sirf_to_torch(tmp_grad, device, requires_grad=True)
        return grad_output*grad, None, None, None

class _AcquisitionModelForward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            sirf_image_template,
            sirf_acq_mdl
            ):
        
        device = torch_image.device
        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        sirf_forward_projected = sirf_acq_mdl.forward(sirf_image_template)
        if torch_image.requires_grad:
            ctx.device = device
            ctx.sirf_forward_projected = sirf_forward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_forward_projected, device, requires_grad=True)
        else:
            return sirf_to_torch(sirf_forward_projected, device)

    @staticmethod
    def backward(ctx,
            grad_output
            ):
        sirf_image = ctx.sirf_acq_mdl.backward(torch_to_sirf(grad_output, ctx.sirf_forward_projected))
        grad = sirf_to_torch(sirf_image, ctx.device, requires_grad=True)
        return grad, None, None, None



class _AcquisitionModelBackward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_measurements,
            sirf_measurements_template,
            sirf_acq_mdl
            ):
        
        device = torch_measurements.device
        sirf_measurements_template = torch_to_sirf(torch_measurements, sirf_measurements_template)
        sirf_backward_projected = sirf_acq_mdl.backward(sirf_measurements_template)
        if torch_measurements.requires_grad:
            ctx.device = device
            ctx.sirf_backward_projected = sirf_backward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_backward_projected, device, requires_grad=True)
        else:
            return sirf_to_torch(sirf_backward_projected, device)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        sirf_measurements = ctx.sirf_acq_mdl.forward(torch_to_sirf(grad_output, ctx.sirf_backward_projected))
        grad = sirf_to_torch(sirf_measurements, ctx.device, requires_grad=True)
        return grad, None, None, None, None

class AcquisitionModelForward(torch.nn.Module):
    def __init__(self,
            acq_mdl, 
            sirf_image_template, 
            sirf_measurements_template, 
            device = "cpu", 
            ):
        super(AcquisitionModelForward, self).__init__()
        # get the shape of image and measurements
        self.acq_mdl = acq_mdl
        self.sirf_image_shape = sirf_image_template.as_array().shape
        self.sirf_measurements_shape = sirf_measurements_template.as_array().shape
        self.sirf_image_template = sirf_image_template
        self.sirf_measurements_template = sirf_measurements_template
        
    def forward(self, image, **kwargs):
        # PyTorch image (2D) is size [batch, channel, height, width] or [batch, height, width]
        # PyTorch volume (3D) is size [batch, channel, depth, height, width] or [batch, depth, height, width]
        # check keys of kwargs
        if not set(kwargs.keys()).issubset({'attenuation', 'background', 'sensitivity','norm','csm'}):
            raise ValueError("Invalid keyword arguments, only 'attenuation', 'background', 'sensitivity','norm','csm' are allowed. Assumes lists of same length.")

        n_batch = image.shape[0]
        if self.sirf_image_shape[0] == 1:
            # if 2D
            if len(image.shape) == 3:
                # if 2D and no channel then add singleton
                # add singleton for SIRF
                image = image.unsqueeze(1)
        else:
            # if 3D
            if len(image.shape) == 4:
                # if 3D and no channel then add singleton
                image = image.unsqueeze(1)
        n_channel = image.shape[1]

        # This looks horrible, but PyTorch should be able to trace.
        batch_images = []
        for batch in n_batch:
            channel_images = []
            for channel in n_channel:
                torch_image = image[batch, channel]
                torch_image = torch_image.view(self.sirf_image_shape)
                # if there are kwargs then raise error
                if kwargs:
                    raise NotImplementedError("Keyword arguments are not implemented yet.")
                channel_images.append(_AcquisitionModelForward.apply(torch_image, self.sirf_image_template, self.sirf_measurements_template, self.acq_mdl))
            batch_images.append(torch.stack(channel_images, dim=0))
        # [batch, channel, *forward_projected.shape]
        return torch.stack(batch_images, dim=0) 
                
if __name__ == '__main__':
    import os
    import numpy
    # Import the PET reconstruction engine
    import sirf.STIR as pet
    # Set the verbosity
    pet.set_verbosity(1)
    # Store temporary sinograms in RAM
    pet.AcquisitionData.set_storage_scheme("memory")
    import sirf
    msg = sirf.STIR.MessageRedirector(info=None, warn=None, errr=None)
    from sirf.Utilities import examples_data_path

    print("2D PET TEST")
    # 2D PET example
    data_path = os.path.join(examples_data_path('PET'), 'thorax_single_slice')
    image = pet.ImageData(os.path.join(data_path,'emission.hv'))
    acq_mdl = pet.AcquisitionModelUsingRayTracingMatrix()
    acq_mdl.set_num_tangential_LORs(5)
    tmp = pet.AcquisitionData(os.path.join(data_path,'template_sinogram.hs'))
    acq_mdl.set_background_term(tmp.get_uniform_copy(0) + 1.5)
    acq_mdl.set_up(tmp,image)
    sirf_measurements = acq_mdl.forward(image)
    # Forward test - check if forward is same between torch and sirf
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    # torch_image,
    torch_image = torch.tensor(image.as_array(), requires_grad=True).to(device)
    # torch_measurements_template,
    torch_measurements_template = torch.tensor(sirf_measurements.as_array(), requires_grad=False).to(device)*0
    # sirf_image_template,
    sirf_image_template = image.get_uniform_copy(0)
    # sirf_acq_mdl
    sirf_acq_mdl = acq_mdl

    
    print("Comparing the forward projected")
    torch_measurements = _AcquisitionModelForward.apply(torch_image, sirf_image_template, sirf_acq_mdl)
    print("Sum of torch: ", torch_measurements.detach().cpu().numpy().sum(), "Sum of sirf: ", sirf_measurements.sum(), "Sum of Differences: ", numpy.abs(torch_measurements.detach().cpu().numpy() - sirf_measurements.as_array()).sum())
    print("Comparing the backward of forward projected")
    # TO TEST THAT WE ARE BACKWARDING CORRECTLY RETAIN GRAD AND SUM THEN BACKWARD
    torch_image.retain_grad()
    torch_measurements.sum().backward()
    # grad sum(Ax) = A^T 1
    comparison = acq_mdl.backward(sirf_measurements.get_uniform_copy(1))
    print("Sum of torch: ", torch_image.grad.sum(), "Sum of sirf: ", comparison.sum(), "Sum of Differences: ", \
        numpy.abs(torch_image.grad.detach().cpu().numpy() - comparison.as_array()).sum())


    torch_measurements = torch.tensor(sirf_measurements.as_array(), requires_grad=True).to(device)
    torch_image_template = torch.tensor(image.as_array(), requires_grad=True).to(device)
    sirf_measurements_template = acq_mdl.forward(image).get_uniform_copy(0)
    bp_sirf_measurements = acq_mdl.backward(sirf_measurements)

    print("Comparing the backward projected")
    torch_image = _AcquisitionModelBackward.apply(torch_measurements, sirf_measurements_template, sirf_acq_mdl)
    print("Sum of torch: ", torch_image.detach().cpu().numpy().sum(), "Sum of sirf: ", bp_sirf_measurements.sum(), \
        "Sum of Differences: ", numpy.abs(torch_image.detach().cpu().numpy() - bp_sirf_measurements.as_array()).sum())
    torch_measurements.retain_grad()
    torch_image.sum().backward()
    # grad sum(A^Tx) = A 1
    comparison = acq_mdl.forward(bp_sirf_measurements.get_uniform_copy(1))
    print("Sum of torch: ", torch_measurements.grad.sum(), "Sum of sirf: ", comparison.sum(), "Sum of Differences: ", \
        numpy.abs(torch_measurements.grad.detach().cpu().numpy() - comparison.as_array()).sum())


    # form objective function
    print("Objective Function Test")
    obj_fun = pet.make_Poisson_loglikelihood(sirf_measurements)
    print("Made poisson loglikelihood")
    obj_fun.set_acquisition_model(sirf_acq_mdl)
    print("Set acquisition model")
    obj_fun.set_up(sirf_image_template)
    print("Set up")

    torch_image = torch.tensor(image.as_array(), requires_grad=True).to(device)
    sirf_image_template = image.get_uniform_copy(0)
    print("Comparing the objective function")
    torch_value = _ObjectiveFunctionModule.apply(torch_image, sirf_image_template, obj_fun) 
    print("Value of torch: ", torch_value.detach().cpu().numpy(), "Value of sirf: ", obj_fun.value(sirf_image_template).as_array())
    torch_image.retain_grad() 
    torch_value.sum().backward()
    # grad sum(f(x)) = f'(x)
    comparison = obj_fun.gradient(sirf_image_template)
    print("Sum of torch: ", torch_image.grad.sum(), "Sum of sirf: ", comparison.sum(), "Sum of Differences: ", (torch_image.grad.detach().cpu().numpy() - comparison.as_array()).sum())
   
    """ from sirf.Gadgetron import *
    print("2D MRI TEST")
    data_path = examples_data_path('MR')
    AcquisitionData.set_storage_scheme('memory')
    input_data = AcquisitionData(
            data_path + '/simulated_MR_2D_cartesian_Grappa2.h5')
    input_norm = input_data.norm()

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    test.check(processed_data.norm() / input_norm, rel_tol=0.01)

    recon = CartesianGRAPPAReconstructor()
    recon.compute_gfactors(False)
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    test.check(complex_images.norm() / input_norm, rel_tol=0.01)

    csms = CoilSensitivityData()

    processed_data.sort()
    csms.calculate(processed_data)

    am = AcquisitionModel(processed_data, complex_images)
    am.set_coil_sensitivity_maps(csms)
    print(am) """

# Objective function
# Acquistion model
#   With torch functionals and comparison
# Jacobian vector for STIR
# similar thing for Gadgetron
# adjointness test
# gradcheck, check traceability