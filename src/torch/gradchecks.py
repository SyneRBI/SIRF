import os
import pytest
import sirf.STIR as pet
import sirf.Gadgetron as mr
from sirf.Utilities import examples_data_path, existing_filepath
import torch
torch.use_deterministic_algorithms(True)
import sirf
from sirf.SIRF_torch import SIRFTorchOperator, SIRFTorchObjectiveFunction, SIRFTorchObjectiveFunctionGradient, sirf_to_torch


def get_data(modality, data_type):
    if modality == "PET":
        pet.set_verbosity(1)
        pet.AcquisitionData.set_storage_scheme("memory")
        msg_pet = pet.MessageRedirector(info=None, warn=None, errr=None)
        pet_data_path = examples_data_path('PET')
        if data_type == "2d":
            raw_data_file = existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
            init_image_file = existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
            acq_data = pet.AcquisitionData(raw_data_file)
            image_data = pet.ImageData(init_image_file)
            acq_model = pet.AcquisitionModelUsingParallelproj()
        elif data_type == "3d":
            raw_data_file = existing_filepath(pet_data_path, 'Utahscat600k_ca_seg4.hs')
            init_image_file = existing_filepath(pet_data_path, 'test_image_PM_QP_6.hv')
            acq_data = pet.AcquisitionData(raw_data_file)
            image_data = pet.ImageData(init_image_file)
            acq_model = pet.AcquisitionModelUsingParallelproj()
        else:
            raise ValueError("Invalid data_type for PET. Choose '2d', '3d'.")
    elif modality == "MR":
        mr.AcquisitionData.set_storage_scheme('memory')
        mr_data_path = examples_data_path('MR')
        if data_type == "2d":
            raw_data_file = existing_filepath(mr_data_path, 'simulated_MR_2D_cartesian.h5')
            input_data = mr.AcquisitionData(raw_data_file)
            prep_gadgets = ['RemoveROOversamplingGadget']
            processed_data = input_data.process(prep_gadgets)
            recon = mr.FullySampledReconstructor()
            recon.set_input(processed_data)
            recon.process()
            complex_images = recon.get_output()
            image_data = complex_images
            csms = mr.CoilSensitivityData()
            processed_data.sort()
            csms.calculate(processed_data)
            acq_model = mr.AcquisitionModel(processed_data, csms)
            acq_model.set_coil_sensitivity_maps(csms)
            acq_data = acq_model.forward(complex_images)
        else:
            raise ValueError("Invalid data_type for MR. Choose '2d'.")
    else:
        raise ValueError("Invalid modality. Only 'PET' or 'MR' are supported.")

    return acq_data, image_data, acq_model


@pytest.fixture(params=[
    ("MR", "2d"), ("PET", "2d"), ("PET", "3d"), 
])
def test_data(request):
    modality, data_type = request.param
    acq_data, image_data, acq_model = get_data(modality, data_type)
    image_data = image_data.get_uniform_copy(1.)
    acq_model.set_up(acq_data, image_data)
    acq_data = acq_model.forward(image_data)
    return acq_data, image_data, acq_model, modality, data_type


test_flags = {
    "forward": True,
    "adjoint": True,
    "objective_wrapped": True,
    "objective": True,
    "objective_gradient": True,
}

def run_gradcheck(func, input_data, data_info, test_name, **kwargs):
    """Helper function to run gradcheck and handle exceptions."""
    modality, data_type = data_info
    try:
        torch.autograd.gradcheck(func, input_data, raise_exception=True, **kwargs)
        print(f"{test_name} passed for {modality} {data_type}")
        assert True
    except Exception as e:
        print(f"{test_name} failed for {modality} {data_type}")
        print(e)
        assert False


@pytest.mark.skipif(not test_flags["forward"], reason="Forward test disabled")
def test_forward_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        pass
    elif modality == "MR":
        pass
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    torch_forward = SIRFTorchOperator(acq_model, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_forward, torch_image_data, (modality, data_type), "Forward",
                  nondet_tol=1e-6, fast_mode=True, eps=1e-2, atol=1e-4, rtol=1e-4)


@pytest.mark.skipif(not test_flags["adjoint"], reason="Adjoint test disabled")
def test_adjoint_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        pass
    elif modality == "MR":
        pass
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    adj_acq_model = sirf.SIRF.AdjointOperator(acq_model)
    torch_adjoint = SIRFTorchOperator(adj_acq_model, acq_data.clone())
    torch_acq_data = torch.tensor(acq_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_adjoint, torch_acq_data, (modality, data_type), "Adjointness",
                  nondet_tol=1e-6, fast_mode=True, eps=1e-2, atol=1e-4, rtol=1e-4)



@pytest.mark.skipif(not test_flags["objective_wrapped"], reason="Wrapped objective test disabled")
def test_objective_function_with_wrapped_acquisition_model_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        # add background term to acq_model
        acq_model.set_background_term(acq_data.get_uniform_copy(1.))
        # add constant term to acq_data 
        acq_data_new = acq_data + 1
        acq_model.set_up(acq_data_new, image_data)
        torch_forward = SIRFTorchOperator(acq_model, image_data.clone())
        torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()
        loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, reduction='sum')
        torch_acq_data = sirf_to_torch(acq_data_new, device=torch_image_data.device)
        torch_acq_model_obj = lambda x: loss(torch_forward(x), torch_acq_data)
    elif modality == "MR":
        torch_forward = SIRFTorchOperator(acq_model, image_data.clone())
        torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()
        torch_acq_data = sirf_to_torch(acq_data, device=torch_image_data.device)
        torch_acq_model_obj = lambda x: (torch_forward(x) - torch_acq_data).abs().pow(2).sum() # L2 squared loss
        
    else:
        pytest.skip("Tests not set up for other modalities at this time.")


    run_gradcheck(torch_acq_model_obj, torch_image_data, (modality, data_type), "Objective (wrapped)",
                  nondet_tol=1e-2, fast_mode=True, eps=1e-3, atol=1e-2, rtol=1e-2)


@pytest.mark.skipif(not test_flags["objective"], reason="Objective test disabled")
def test_objective_function_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        # add background term to acq_model
        acq_model.set_background_term(acq_data.get_uniform_copy(1.))
        # add constant term to acq_data 
        acq_data_new = acq_data + 1
        obj_fun = pet.make_Poisson_loglikelihood(acq_data_new)
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_up(image_data)

    torch_obj_fun = SIRFTorchObjectiveFunction(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()


    run_gradcheck(torch_obj_fun, torch_image_data, (modality, data_type), "Objective",
                  nondet_tol=1e-2, fast_mode=True, eps=1e-3, atol=1e-2, rtol=1e-2)


@pytest.mark.skipif(not test_flags["objective_gradient"], reason="Objective Gradient test disabled")
def test_objective_function_gradient_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        # add background term to acq_model
        acq_model.set_background_term(acq_data.get_uniform_copy(1.))
        # add constant term to acq_data 
        acq_data_new = acq_data + 1
        obj_fun = pet.make_Poisson_loglikelihood(acq_data_new)
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_up(image_data)


    torch_obj_fun_grad = SIRFTorchObjectiveFunctionGradient(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_obj_fun_grad, torch_image_data, (modality, data_type), "Objective Gradient",
                  nondet_tol=1e-4, fast_mode=True, eps=1e-3, atol=1e-2, rtol=1e-2)