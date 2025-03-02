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
import matplotlib.pyplot as plt
import torch
# set deterministic computing
torch.use_deterministic_algorithms(True)

from SIRF_torch import *
# -	Gradient checks
# -	Adjointness checks

def get_pet_2d():
    pet_data_path = examples_data_path('PET')
    pet_2d_raw_data_file = existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
    pet_2d_acq_data = AcquisitionData(pet_2d_raw_data_file)
    pet_2d_init_image_file = existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
    pet_2d_image_data = ImageData(pet_2d_init_image_file)
    pet_2d_acq_model = AcquisitionModelUsingParallelproj()
    """ pet_2d_acq_model = AcquisitionModelUsingRayTracingMatrix()
    pet_2d_acq_model.set_num_tangential_LORs(5) """
    pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)
    pet_2d_acq_model.adjoint(pet_2d_acq_data)
    return pet_2d_acq_data, pet_2d_image_data, pet_2d_acq_model

def get_pet_3d():
    pet_data_path = examples_data_path('PET')
    pet_3d_raw_data_file = existing_filepath(pet_data_path, 'Utahscat600k_ca_seg4.hs')
    pet_3d_acq_data = AcquisitionData(pet_3d_raw_data_file)
    pet_3d_init_image_file = existing_filepath(pet_data_path, 'test_image_PM_QP_6.hv')
    pet_3d_image_data = ImageData(pet_3d_init_image_file)
    pet_3d_acq_model = AcquisitionModelUsingParallelproj()
    """ pet_3d_acq_model = AcquisitionModelUsingRayTracingMatrix()
    pet_3d_acq_model.set_num_tangential_LORs(5) """
    pet_3d_acq_model.set_background_term(pet_2d_acq_data.get_uniform_copy(0) + 1.5)
    return pet_3d_acq_data, pet_3d_image_data, pet_3d_acq_model

def rescale_image(image_data, scale = 8):
    new_size = (image_data.shape[0]//scale, image_data.shape[1]//scale, image_data.shape[2]//scale)
    image_data = image_data.zoom_image(zooms=(1/scale, 1/scale, 1/scale), offsets_in_mm=(0.0, 0.0, 0.0), size=new_size)
    return image_data
    

if __name__ == "__main__":
    # https://pytorch.org/docs/stable/notes/gradcheck.html
    acq_data, image_data, acq_model = get_pet_2d()
    acq_data = acq_data.get_uniform_copy(1.)
    image_data = image_data.get_uniform_copy(1.)
    acq_model.set_up(acq_data, image_data)

    torch_forward = SIRFOperator(acq_model, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad = True).unsqueeze(0).cuda()

    try:
        torch.autograd.gradcheck(torch_forward,
            torch_image_data,
            nondet_tol=1e-6, 
            fast_mode=True, 
            eps=1e-3,
            atol=1e-4,
            rtol=1e-4,
            raise_exception=True
            )
        print("Forward gradcheck passed")
    except:
        print("Forward gradcheck failed")

    adj_acq_model = sirf.SIRF.AdjointOperator(acq_model)
    torch_adjoint = SIRFOperator(adj_acq_model, acq_data.clone())
    torch_acq_data = torch.tensor(acq_data.as_array(), requires_grad = True).unsqueeze(0).cuda()
    
    try:
        torch.autograd.gradcheck(torch_adjoint,
            torch_acq_data,
            nondet_tol=1e-6, 
            fast_mode=True, 
            eps=1e-3,
            atol=1e-4,
            rtol=1e-4,
            raise_exception=True
            )
        print("Adjointness gradcheck passed")
    except:
        print("Adjointness gradcheck failed")

    acq_model.set_background_term(acq_data.get_uniform_copy(1.))
    acq_model.set_up(acq_data, image_data)

    loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum')
    torch_acq_data = sirf_to_torch(acq_data, device = torch_image_data.device)  
    torch_acq_model_obj = lambda x: loss(torch_forward(x), torch_acq_data)

    try:
        torch.autograd.gradcheck(torch_acq_model_obj,
            torch_image_data,
            nondet_tol=1e-4, 
            fast_mode=True, 
            eps=1e-3,
            atol=1e-2,
            rtol=1e-2,
            raise_exception=True
            )
        print("Objective function with wrapped aquisition model gradcheck passed")
    except:
        print("Objective function with wrapped aquisition model gradcheck failed")

    # set up the objective function
    obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_up(image_data)

    # test the objective function
    torch_obj_fun = SIRFObjectiveFunction(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad = True).unsqueeze(0).cuda()

    try:
        torch.autograd.gradcheck(torch_obj_fun,
            torch_image_data,
            nondet_tol=1e-4, 
            fast_mode=True, 
            eps=1e-3,
            atol=1e-2,
            rtol=1e-2,
            raise_exception=True,
            )
        print("Objective function gradcheck passed")
    except:
        print("Objective function gradcheck failed")

    
    # test the objective function
    torch_obj_fun_grad = SIRFObjectiveFunctionGradient(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad = True).unsqueeze(0).cuda()

    
    try:
        torch.autograd.gradcheck(torch_obj_fun_grad,
            torch_image_data,
            nondet_tol=1e-6, 
            fast_mode=True, 
            eps=1e-3,
            atol=1e-4,
            rtol=1e-4,
            raise_exception=True,
                )
        print("Objective function gradient gradcheck passed")
    except:
        print("Objective function gradient gradcheck failed")