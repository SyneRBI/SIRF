
from SIRF_torch import *
import sirf.STIR as pet
pet.set_verbosity(1)
pet.AcquisitionData.set_storage_scheme("memory")
import sirf
msg = sirf.STIR.MessageRedirector(info=None, warn=None, errr=None)
from sirf.Utilities import examples_data_path
import matplotlib.pyplot as plt

def get_pet_2d():
    pet_data_path = examples_data_path('PET')
    pet_2d_raw_data_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
    pet_2d_acq_data = pet.AcquisitionData(pet_2d_raw_data_file)
    pet_2d_init_image_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
    pet_2d_image_data = pet.ImageData(pet_2d_init_image_file)
    pet_2d_acq_model = pet.AcquisitionModelUsingParallelproj()
    pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)
    pet_2d_acq_data = pet_2d_acq_model.forward(pet_2d_image_data) + 1.5
    pet_2d_acq_model.set_background_term(pet_2d_acq_data.get_uniform_copy(0) + 1.5)
    pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)
    return pet_2d_acq_data, pet_2d_image_data, pet_2d_acq_model

def get_pet_3d():
    pet_data_path = examples_data_path('PET')
    pet_3d_raw_data_file = pet.existing_filepath(pet_data_path, 'Utahscat600k_ca_seg4.hs')
    pet_3d_acq_data = pet.AcquisitionData(pet_3d_raw_data_file)
    pet_3d_init_image_file = pet.existing_filepath(pet_data_path, 'test_image_PM_QP_6.hv')
    pet_3d_image_data = pet.ImageData(pet_3d_init_image_file)
    pet_3d_acq_model = pet.AcquisitionModelUsingRayTracingMatrix()
    pet_3d_acq_model.set_num_tangential_LORs(5)
    pet_2d_acq_model = pet.AcquisitionModelUsingParallelproj()
    pet_3d_acq_model.set_up(pet_3d_acq_data, pet_3d_image_data)
    return pet_3d_acq_data, pet_3d_image_data, pet_3d_acq_model


# LEARNED PRIMAL DUAL


# UNROLLED VARNET-like network - backpropagation through gradient of objective


# Gradient Descent with ADAM (wrapped objective and acq model)


class ConvBlock(torch.nn.Module):
    def __init__(self, in_channels=2, out_channels=1, kernel_size=3, padding=1, hidden_channels=10):
        super(ConvBlock, self).__init__()
        self.conv1 = torch.nn.Conv2d(in_channels, hidden_channels, kernel_size, padding=padding)
        self.conv2 = torch.nn.Conv2d(hidden_channels, hidden_channels, kernel_size, padding=padding)
        self.conv3 = torch.nn.Conv2d(hidden_channels, out_channels, kernel_size, padding=padding)
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = self.relu(self.conv2(x))
        x = self.conv3(x)
        return x


class PETVarNet(torch.nn.Module):
    def __init__(self, ConvBlocks, ObjectiveFunctionGradient):
        super(PETVarNet, self).__init__()
        self.ConvBlocks = torch.nn.ModuleList(ConvBlocks)
        self.ObjectiveFunctionGradient = ObjectiveFunctionGradient
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        for i in range(len(self.ConvBlocks)):
            filtered_x = self.relu(self.ConvBlocks[i](x))
            x = filtered_x + self.relu(self.ObjectiveFunctionGradient(filtered_x))
        return x

class PETLearnedPrimalDual(torch.nn.Module):
    def __init__(self, PrimalConvBlocks, DualConvBlocks, FwdOperator, AdjOperator):
        super(PETLearnedPrimalDual, self).__init__()
        # check the length of primal and dual conv blocks
        if len(PrimalConvBlocks) != len(DualConvBlocks):
            raise ValueError("Length of Primal and Dual Conv Blocks must be the same")
        self.PrimalConvBlocks = torch.nn.ModuleList(PrimalConvBlocks)
        self.DualConvBlocks = torch.nn.ModuleList(DualConvBlocks)
        # add sirf singleton for channel then remove sirf singletons acq data
        self.FwdOperator = lambda x: FwdOperator(x).squeeze(2).squeeze(3)
        # add sirf singleton for channel then remove sirf singleton 2D image
        self.AdjOperator = lambda x: AdjOperator(x.unsqueeze(2).unsqueeze(2)).squeeze(2)
        self.relu = torch.nn.ReLU()
    
    def forward(self, y):
        batch_size = y.shape[0]
        h = torch.zeros_like(y, device=y.device)
        dummy_f = self.AdjOperator(y)
        f = torch.zeros_like(dummy_f, device=y.device)
        n_iter = len(self.PrimalConvBlocks)

        for i in range(n_iter):
            Op_f = self.FwdOperator(f)
            dual_input = torch.cat([Op_f, y], dim=1)
            h = self.DualConvBlocks[i](dual_input)
            h = self.relu(h)
            OpAdj_h = self.AdjOperator(h)
            f = self.PrimalConvBlocks[i](torch.cat([f,OpAdj_h], dim=1)) 
            f = self.relu(f) 
        return f




class UseCases:
    def __init__(self, modality, dim, test_name):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        if modality == 'PET':
            if dim == '2D':
                self.acq_data, self.image, self.acq_model = get_pet_2d()
                # make torch versions of acq_data and image
                self.torch_acq_data = sirf_to_torch(self.acq_data, self.device)
                self.torch_image = sirf_to_torch(self.image, self.device).unsqueeze(0)
            elif dim == '3D':
                self.acq_data, self.image, self.acq_model = get_pet_3d()
        if test_name == 'all':
            if dim == '2D':
                self.pet_varnet()
                self.learned_primal_dual()
            self.compare_gradient_descents()
        elif test_name == 'gradient_descent':
            self.compare_gradient_descents()
        elif test_name == 'unrolled_backward':
            self.unrolled_backward_test()
        elif test_name == 'objective_function':
            self.objective_function_test()



    def learned_primal_dual(self):
        # use seed for reproducibility

        # set up learned primal dual network
        PrimalConvBlocks = [ConvBlock() for i in range(3)]
        DualConvBlocks = [ConvBlock() for i in range(3)]
        # set up the forward and adjoint operators
        sirf_image_template = self.image.get_uniform_copy(1)
        FwdOperator = SIRFTorchOperator(self.acq_model, self.image.get_uniform_copy(1))
        AdjOperator = SIRFTorchOperator(sirf.SIRF.AdjointOperator(self.acq_model.get_linear_acquisition_model()), self.acq_data.get_uniform_copy(1))
        torch_target = sirf_to_torch(sirf_image_template, self.device).unsqueeze(0)
        net = PETLearnedPrimalDual(PrimalConvBlocks, DualConvBlocks, FwdOperator, AdjOperator)
        net.to(self.device)

        # build the objective function
        relu = torch.nn.ReLU()
        loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum')
        torch_input = sirf_to_torch(self.acq_data, self.device)
        torch_obj_func = lambda x: -loss(FwdOperator(relu(net(x))), torch_input)

        # set up the optimizer
        optimizer = torch.optim.Adam(net.parameters(), lr=2e-3)
        for i in range(5):
            optimizer.zero_grad()
            loss_val = torch_obj_func(torch_input)
            loss_val.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss_val.item())
        out = net(torch_input)
        plt.imshow(out.detach().cpu().numpy()[0,0])
        plt.savefig("learned_primal_dual.png")
        

    def pet_varnet(self):
        # use seed for reproducibility
        torch.manual_seed(42)

        # set up pet varnet network
        ConvBlocks = [ConvBlock(in_channels=1) for i in range(3)]
        # set up the Objective Function Gradient
        obj_fun = pet.make_Poisson_loglikelihood(self.acq_data)
        print("Made poisson loglikelihood")
        obj_fun.set_acquisition_model(self.acq_model)
        print("Set acquisition model")
        obj_fun.set_up(self.image.get_uniform_copy(1))
        print("Set up")
        ObjectiveFunctionGradient = SIRFTorchObjectiveFunctionGradient(obj_fun, self.image.get_uniform_copy(1))
        print("Made objective function gradient")
        net = PETVarNet(ConvBlocks, ObjectiveFunctionGradient)
        net.to(self.device)


        # build the objective function
        ObjectiveFunction = SIRFTorchObjectiveFunction(obj_fun, self.image.get_uniform_copy(1))
        # non-negative image as input
        torch_input = torch.ones_like(self.torch_image, device=self.device)
        relu = torch.nn.ReLU()

        # set up the optimizer
        optimizer = torch.optim.Adam(net.parameters(), lr=2e-3)
        for i in range(20):
            optimizer.zero_grad()
            loss_val = - ObjectiveFunction(relu(net(torch_input)))
            loss_val.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss_val.item())
        out = relu(net(torch_input))
        plt.imshow(out.detach().cpu().numpy()[0,0])
        plt.savefig("pet_varnet.png")


    def compare_gradient_descents(self):
        print("Comparing Gradient Descents")
        lr = 0.1
        n_iter = 3
        print(f"Learning Rate: {lr}, Number of Iterations: {n_iter}")
        print("Gradient Descent with Acquisition Model")
        out_acq_model = self.gradient_descent_with_acq_model(lr=lr, n_iter=n_iter)
        print("Gradient Descent with Objective Function")
        out_obj_func = self.gradient_descent_with_obj_func(lr=lr, n_iter=n_iter)
        print("Sum of Differences: ", numpy.abs(out_acq_model - out_obj_func).sum())

    def gradient_descent_with_acq_model(self, lr, n_iter):
        # PARAMETERS ARE IMAGE-LIKE
        torch_image = sirf_to_torch(self.image.get_uniform_copy(1), self.device).unsqueeze(0) # add bactch dimension
        torch_image_init = torch_image.clone()
        torch_image_params = torch.nn.Parameter(torch_image)
        torch_measurements = sirf_to_torch(self.acq_data, self.device)
        torch_acq_model = SIRFTorchOperator(self.acq_model, self.image.get_uniform_copy(1))
        relu = torch.nn.ReLU()
        loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum')
        torch_obj_func = lambda x: loss(torch_acq_model(relu(x)), torch_measurements)
        optimizer = torch.optim.Adam([torch_image_params], lr=lr)
        for i in range(n_iter):
            optimizer.zero_grad()
            loss_val = torch_obj_func(torch_image_params)
            loss_val.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss_val.item())
        return torch_image_params.data.detach().cpu().numpy()

    def gradient_descent_with_obj_func(self, lr, n_iter):
        # PARAMETERS ARE IMAGE-LIKE
        torch_image = sirf_to_torch(self.image.get_uniform_copy(1), self.device).unsqueeze(0)
        torch_image_init = torch_image.clone()
        torch_image_params = torch.nn.Parameter(torch_image)
        obj_fun = pet.make_Poisson_loglikelihood(self.acq_data)
        print("Made poisson loglikelihood")
        obj_fun.set_acquisition_model(self.acq_model)
        print("Set acquisition model")
        obj_fun.set_up(self.image.get_uniform_copy(1))
        print("Set up")
        obj = SIRFTorchObjectiveFunction(obj_fun, self.image.get_uniform_copy(1))
        relu = torch.nn.ReLU()
        obj_relu = lambda x: obj(relu(x))
        optimizer = torch.optim.Adam([torch_image_params], lr=lr)
        for i in range(n_iter):
            optimizer.zero_grad()
            loss = obj(torch_image_params)
            loss.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss.item())
        return torch_image_params.data.detach().cpu().numpy()
        


if __name__ == '__main__':
    torch.manual_seed(42)
    uses = UseCases('PET', '2D', 'all')