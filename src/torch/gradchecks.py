import os
import pytest
import sirf.STIR as pet
from sirf.Utilities import examples_data_path, existing_filepath
import torch
from SIRF_torch import *
import sirf
# Set verbosity and storage scheme
pet.set_verbosity(1)
pet.AcquisitionData.set_storage_scheme("memory")
torch.use_deterministic_algorithms(True)


def get_data(modality, data_type):
    if modality == "PET":
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
    else:
        raise ValueError("Invalid modality. Only 'PET' at the moment.")

    return acq_data, image_data, acq_model


@pytest.fixture(params=[
    ("PET", "2d"), ("PET", "3d")
])
def test_data(request):
    modality, data_type = request.param
    acq_data, image_data, acq_model = get_data(modality, data_type)
    acq_data = acq_data.get_uniform_copy(1.)
    image_data = image_data.get_uniform_copy(1.)
    acq_model.set_up(acq_data, image_data)
    return acq_data, image_data, acq_model, modality, data_type


test_flags = {
    "forward": False,
    "adjoint": False,
    "objective_wrapped": False,
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
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    torch_forward = SIRFTorchOperator(acq_model, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_forward, torch_image_data, (modality, data_type), "Forward",
                  nondet_tol=1e-6, fast_mode=True, eps=1e-3, atol=1e-4, rtol=1e-4)


@pytest.mark.skipif(not test_flags["adjoint"], reason="Adjoint test disabled")
def test_adjoint_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        pass
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    adj_acq_model = sirf.SIRF.AdjointOperator(acq_model)
    torch_adjoint = SIRFTorchOperator(adj_acq_model, acq_data.clone())
    torch_acq_data = torch.tensor(acq_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_adjoint, torch_acq_data, (modality, data_type), "Adjointness",
                  nondet_tol=1e-6, fast_mode=True, eps=1e-3, atol=1e-4, rtol=1e-4)



@pytest.mark.skipif(not test_flags["objective_wrapped"], reason="Wrapped objective test disabled")
def test_objective_function_with_wrapped_acquisition_model_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        pass
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    acq_model.set_background_term(acq_data.get_uniform_copy(1.))
    acq_model.set_up(acq_data, image_data)

    torch_forward = SIRFTorchOperator(acq_model, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()
    loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, reduction='sum')
    torch_acq_data = sirf_to_torch(acq_data, device=torch_image_data.device)
    torch_acq_model_obj = lambda x: loss(torch_forward(x), torch_acq_data)

    run_gradcheck(torch_acq_model_obj, torch_image_data, (modality, data_type), "Objective (wrapped)",
                  nondet_tol=1e-2, fast_mode=True, eps=1e-3, atol=1e-2, rtol=1e-2)


@pytest.mark.skipif(not test_flags["objective"], reason="Objective test disabled")
def test_objective_function_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_up(image_data)

    torch_obj_fun = SIRFTorchObjectiveFunction(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()


    run_gradcheck(torch_obj_fun, torch_image_data, (modality, data_type), "Objective",
                  nondet_tol=1e-1, fast_mode=True, eps=1e-3, atol=1e-1, rtol=1e-1)


@pytest.mark.skipif(not test_flags["objective_gradient"], reason="Objective Gradient test disabled")
def test_objective_function_gradient_gradcheck(test_data):
    acq_data, image_data, acq_model, modality, data_type = test_data
    if modality == "PET":
        obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    else:
        pytest.skip("Tests not set up for other modalities at this time.")

    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_up(image_data)


    torch_obj_fun_grad = SIRFTorchObjectiveFunctionGradient(obj_fun, image_data.clone())
    torch_image_data = torch.tensor(image_data.as_array(), requires_grad=True).unsqueeze(0).cuda()

    run_gradcheck(torch_obj_fun_grad, torch_image_data, (modality, data_type), "Objective Gradient",
                  nondet_tol=1e-6, fast_mode=True, eps=1e-3, atol=1e-4, rtol=1e-4)