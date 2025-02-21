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


from SIRF_torch import *
# -	LPD
# -	DIP
# -	ADAM optimisation of SIRF objective function
# - PNLL
# -	Torch OSEM implementation

def get_pet_2d():
    pet_data_path = examples_data_path('PET')
    pet_2d_raw_data_file = existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
    pet_2d_acq_data = AcquisitionData(pet_2d_raw_data_file)
    pet_2d_init_image_file = existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
    pet_2d_image_data = ImageData(pet_2d_init_image_file)
    pet_2d_acq_mdl = AcquisitionModelUsingRayTracingMatrix()
    pet_2d_acq_mdl.set_num_tangential_LORs(5)
    pet_2d_acq_mdl.set_background_term(pet_2d_acq_data.get_uniform_copy(0) + 1.5)
    pet_2d_acq_mdl.set_up(pet_2d_acq_data, pet_2d_image_data)
    pet_2d_acq_data = pet_2d_acq_mdl.forward(pet_2d_image_data)
    return pet_2d_acq_data, pet_2d_image_data, pet_2d_acq_mdl

def get_pet_3d():
    pet_data_path = examples_data_path('PET')
    pet_3d_raw_data_file = existing_filepath(pet_data_path, 'Utahscat600k_ca_seg4.hs')
    pet_3d_acq_data = AcquisitionData(pet_3d_raw_data_file)
    pet_3d_init_image_file = existing_filepath(pet_data_path, 'test_image_PM_QP_6.hv')
    pet_3d_image_data = ImageData(pet_3d_init_image_file)
    pet_3d_acq_mdl = AcquisitionModelUsingRayTracingMatrix()
    pet_3d_acq_mdl.set_num_tangential_LORs(5)
    pet_3d_acq_mdl.set_up(pet_3d_acq_data, pet_3d_image_data)
    return pet_3d_acq_data, pet_3d_image_data, pet_3d_acq_mdl


class unrolled_fwd(torch.nn.Module):
    def __init__(
            self, 
            acq_mdl, 
            sirf_image_template, 
            sirf_measurements_template, 
            device = "cpu"):
        super(unrolled_fwd, self).__init__()
        self.wrapped_forward = SIRFAcquisitionModelForward(acq_mdl, sirf_image_template, sirf_measurements_template, device)
        self.ReLU = torch.nn.ReLU()

    def forward(self, x):
        return self.wrapped_forward(self.ReLU(x))

class unrolled_bwd(torch.nn.Module):
    def __init__(
            self, 
            acq_mdl, 
            sirf_image_template, 
            sirf_measurements_template, 
            device = "cpu"):
        super(unrolled_bwd, self).__init__()
        self.wrapped_backward = SIRFAcquisitionModelBackward(acq_mdl, sirf_image_template, sirf_measurements_template, device)
        self.ReLU = torch.nn.ReLU()


    def forward(self, x):
        return self.wrapped_backward(self.ReLU(x))


class test_wrapper:
    def __init__(self, modality, dim, test_name):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        if modality == 'PET':
            if dim == '2D':
                self.acq_data, self.image, self.acq_mdl = get_pet_2d()
            elif dim == '3D':
                self.acq_data, self.image, self.acq_mdl = get_pet_3d()
        if test_name == 'all':
            self.unrolled_forward_test()
            #self.unrolled_backward_test()
            self.objective_function_test()
        elif test_name == 'unrolled_forward':
            self.unrolled_forward_test()
        elif test_name == 'unrolled_backward':
            self.unrolled_backward_test()
        elif test_name == 'objective_function':
            self.objective_function_test()

    def unrolled_forward_test(self):
        # PARAMETERS ARE IMAGE-LIKE
        print("Unrolled Forward Test")
        torch_image = sirf_to_torch(self.image.get_uniform_copy(1), self.device).unsqueeze(0) # add bactch dimension
        torch_image_init = torch_image.clone()
        torch_image_params = torch.nn.Parameter(torch_image)
        model = unrolled_fwd(self.acq_mdl, self.image.get_uniform_copy(1), self.acq_data, self.device)
        # make lambda function for objective function
        data = sirf_to_torch(self.acq_data, self.device)
        loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum')
        obj = lambda x: loss(model(x), data)
        self.optimise_with_adam([torch_image_params], obj, lr=0.1, n_iter=100)
        # plot before and after, and target
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        fig.colorbar(axs[0].imshow(torch_image_init[0, 0, :, :].cpu().numpy()), ax=axs[0])
        axs[0].set_title("Initial Image")
        fig.colorbar(axs[1].imshow(torch_image_params.data.detach().cpu().numpy()[0, 0, :, :]), ax=axs[1])
        axs[1].set_title("Reconstructed Image")
        fig.colorbar(axs[2].imshow(self.image.as_array()[0, :, :]), ax=axs[2])
        axs[2].set_title("Target Image")
        plt.savefig("unrolled_forward_test.png")

    def unrolled_backward_test(self):
        # PARAMETERS ARE MEASUREMENT-LIKE
        print("Unrolled Backward Test")
        torch_measurements = sirf_to_torch(self.acq_data, self.device)
        torch_measurements_params = torch.nn.Parameter(torch_measurements)
        net = unrolled_bwd(self.acq_mdl, self.image, self.acq_data, self.device)
        # make lambda function for objective function
        image = sirf_to_torch(self.image, self.device)
        obj = lambda x: (net(x) - image).pow(2).sum()
        self.optimise_with_adam([torch_measurements_params], obj)

    def objective_function_test(self):
        # PARAMETERS ARE IMAGE-LIKE
        print("Objective Function Test")
        torch_image = sirf_to_torch(self.image.get_uniform_copy(0), self.device).unsqueeze(0)
        torch_image_init = torch_image.clone()
        torch_image_params = torch.nn.Parameter(torch_image)
        obj_fun = pet.make_Poisson_loglikelihood(self.acq_data)
        obj_fun.set_acquisition_model(self.acq_mdl)
        obj_fun.set_up(self.image)
        obj = SIRFObjectiveFunction(obj_fun, self.image.get_uniform_copy(1), self.device)
        relu = torch.nn.ReLU()
        obj_relu = lambda x: obj(relu(x))
        self.optimise_with_adam([torch_image_params], obj, lr=.5, n_iter=100)
        # plot before and after, and target
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        fig.colorbar(axs[0].imshow(torch_image_init[0, 0, :, :].cpu().numpy()), ax=axs[0])
        axs[0].set_title("Initial Image")
        fig.colorbar(axs[1].imshow(torch_image_params.data.detach().cpu().numpy()[0, 0, :, :]), ax=axs[1])
        axs[1].set_title("Reconstructed Image")
        fig.colorbar(axs[2].imshow(self.image.as_array()[0, :, :]), ax=axs[2])
        axs[2].set_title("Target Image")
        plt.savefig("objective_function_test.png")
        

    def optimise_with_adam(self, params, obj, lr=0.01, n_iter=100):
        optimizer = torch.optim.Adam(params, lr=lr)
        for i in range(n_iter):
            optimizer.zero_grad()
            loss = obj(params[0])
            loss.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss.item())
            if i % 10 == 0:
                print("Iteration: ", i, "Loss: ", loss.item())

test = test_wrapper('PET', '2D', 'all')