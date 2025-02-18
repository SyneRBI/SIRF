try:
    import torch
except ModuleNotFoundError:
    raise ModuleNotFoundError('Failed to import torch. Please install PyTorch first.')


# based on 
# https://github.com/educating-dip/pet_deep_image_prior/blob/main/src/deep_image_prior/torch_wrapper.py


def sirf_to_torch(sirf_src, torch_dest, requires_grad=False):
    if requires_grad:
        return torch.tensor(sirf_src.as_array(), requires_grad=True).to(torch_dest.device)
    else:
        return torch.tensor(sirf_src.as_array()).to(torch_dest.device)

def torch_to_sirf(torch_src, sirf_dest): 
    return sirf_dest.fill(torch_src.detach().cpu().numpy())


class _ObjectiveFunctionModule(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            sirf_image_template, 
            sirf_obj_func
            ):

        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        value_np = sirf_obj_func.get_value(sirf_image_template).as_array()
        if torch_image.requires_grad:
            ctx.save_for_backward(torch_image, sirf_image_template, sirf_obj_func)
            return torch.tensor(value_np, requires_grad=True).to(torch_image.device)
        else:
            return torch.tensor(value_np).to(torch_image.device)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        torch_image, sirf_image_template, sirf_obj_func = ctx.saved_tensors
        tmp_grad = sirf_obj.get_gradient(sirf_image_template)
        grad = sirf_to_torch(tmp_grad, torch_image, requires_grad=True)
        return grad_output*grad, None, None, None

class _AcquisitionModelForward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            torch_measurements_template,
            sirf_image_template,
            sirf_acq_mdl
            ):
        
        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        sirf_forward_projected = sirf_acq_mdl.forward(sirf_image_template)
        if torch_image.requires_grad:
            ctx.torch_image = torch_image
            ctx.sirf_forward_projected = sirf_forward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_forward_projected, torch_measurements_template, requires_grad=True)
        else:
            return sirf_to_torch(sirf_forward_projected, torch_measurements_template)

    @staticmethod
    def backward(ctx,
            grad_output
            ):
        sirf_image = ctx.sirf_acq_mdl.backward(torch_to_sirf(grad_output, ctx.sirf_forward_projected))
        grad = sirf_to_torch(sirf_image, ctx.torch_image, requires_grad=True)
        return grad, None, None, None, None



class _AcquisitionModelBackward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_measurements,
            torch_image_template,
            sirf_measurements_template,
            sirf_acq_mdl
            ):
        
        sirf_measurements_template = torch_to_sirf(torch_measurements, sirf_measurements_template)
        sirf_backward_projected = sirf_acq_mdl.backward(sirf_measurements_template)
        if torch_image_template.requires_grad:
            ctx.torch_measurements = torch_measurements
            ctx.sirf_backward_projected = sirf_backward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_backward_projected, torch_image_template, requires_grad=True)
        else:
            return sirf_to_torch(sirf_backward_projected, torch_image_template)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        sirf_measurements = ctx.sirf_acq_mdl.forward(torch_to_sirf(grad_output, ctx.sirf_backward_projected))
        grad = sirf_to_torch(sirf_measurements, ctx.torch_measurements, requires_grad=True)
        return grad, None, None, None, None


if __name__ == '__main__':
    import os
    import numpy
    import sirf.STIR as pet
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
    torch_measurements = _AcquisitionModelForward.apply(torch_image, torch_measurements_template, sirf_image_template, sirf_acq_mdl)
    print("Sum of torch: ", torch_measurements.detach().cpu().numpy().sum(), "Sum of sirf: ", sirf_measurements.sum(), "Sum of Differences: ", (torch_measurements.detach().cpu().numpy() - sirf_measurements.as_array()).sum())
    print("Comparing the backward of forward projected")
    # TO TEST THAT WE ARE BACKWARDING CORRECTLY RETAIN GRAD AND SUM THEN BACKWARD
    torch_image.retain_grad()
    torch_measurements.sum().backward()
    # grad sum(Ax) = A^T 1
    comparison = acq_mdl.backward(sirf_measurements.get_uniform_copy(1))
    print("Sum of torch: ", torch_image.grad.sum(), "Sum of sirf: ", comparison.sum(), "Sum of Differences: ", (torch_image.grad.detach().cpu().numpy() - comparison.as_array()).sum())


    torch_measurements = torch.tensor(sirf_measurements.as_array(), requires_grad=True).to(device)
    torch_image_template = torch.tensor(image.as_array(), requires_grad=True).to(device)
    sirf_measurements_template = acq_mdl.forward(image).get_uniform_copy(0)
    bp_sirf_measurements = acq_mdl.backward(sirf_measurements)

    print("Comparing the backward projected")
    torch_image = _AcquisitionModelBackward.apply(torch_measurements, torch_image_template, sirf_measurements_template, sirf_acq_mdl)
    print("Sum of torch: ", torch_image.detach().cpu().numpy().sum(), "Sum of sirf: ", bp_sirf_measurements.sum(), "Sum of Differences: ", (torch_image.detach().cpu().numpy() - bp_sirf_measurements.as_array()).sum())
    torch_measurements.retain_grad()
    torch_image.sum().backward()
    # grad sum(A^Tx) = A 1
    comparison = acq_mdl.forward(bp_sirf_measurements.get_uniform_copy(1))
    print("Sum of torch: ", torch_measurements.grad.sum(), "Sum of sirf: ", comparison.sum(), "Sum of Differences: ", (torch_measurements.grad.detach().cpu().numpy() - comparison.as_array()).sum())


    # form objective function
    print("Objective Function Test")
    obj_fun = pet.make_Poisson_loglikelihood(sirf_measurements)
    obj_fun.set_acquisition_model(acq_mdl)
    obj_fun.set_up(sirf_image_template)

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



# Objective function
# Acquistion model
#   With torch functionals and comparison
# Jacobian vector for STIR
# similar thing for Gadgetron
# adjointness test
# gradcheck, check traceability