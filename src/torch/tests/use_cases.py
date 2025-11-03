from sirf.torch import Operator, ObjectiveFunction, ObjectiveFunctionGradient, sirf_to_torch
import sirf.STIR as pet
from sirf.Utilities import examples_data_path
import sirf.SIRF as sirf
import torch
import numpy as np
import matplotlib.pyplot as plt

# Set global settings for STIR
pet.set_verbosity(1)
pet.AcquisitionData.set_storage_scheme("memory")

# Redirect STIR messages to avoid cluttering the output (optional)
msg = pet.MessageRedirector(info=None, warn=None, errr=None)


def get_pet_2d():
    """
    Loads 2D PET data, sets up the acquisition model, and adds a background term.

    Returns:
        tuple: A tuple containing the acquisition data, initial image, and the acquisition model.
    """
    pet_data_path = examples_data_path('PET')
    pet_2d_raw_data_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
    pet_2d_acq_data = pet.AcquisitionData(pet_2d_raw_data_file)

    pet_2d_init_image_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
    pet_2d_image_data = pet.ImageData(pet_2d_init_image_file)

    pet_2d_acq_model = pet.AcquisitionModelUsingParallelproj()
    pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)

    # Simulate acquisition data and add a background term (for realism)
    pet_2d_acq_data = pet_2d_acq_model.forward(pet_2d_image_data) + 5.
    pet_2d_acq_model.set_background_term(pet_2d_acq_data.get_uniform_copy(5.))

    # Truncate the image to a cylinder (optional, for focusing on a region)
    data_processor = pet.TruncateToCylinderProcessor()
    pet_2d_acq_model.set_image_data_processor(data_processor)
    pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)

    return pet_2d_acq_data, pet_2d_image_data, pet_2d_acq_model


class PETLearnedPrimalDual(torch.nn.Module):
    """
    A PyTorch module implementing a learned primal-dual algorithm for PET image reconstruction.
    """
    def __init__(self, PrimalConvBlocks, DualConvBlocks, FwdOperator, AdjOperator):
        """
        Initializes the PETLearnedPrimalDual module.

        Args:
            PrimalConvBlocks (list): A list of convolutional blocks for the primal update.
            DualConvBlocks (list): A list of convolutional blocks for the dual update.
            FwdOperator (callable): A callable object representing the forward operator (A).
            AdjOperator (callable): A callable object representing the adjoint operator (A^T).
        """
        super(PETLearnedPrimalDual, self).__init__()

        # Check if the number of primal and dual blocks is the same
        if len(PrimalConvBlocks) != len(DualConvBlocks):
            raise ValueError("Length of Primal and Dual Conv Blocks must be the same")

        self.PrimalConvBlocks = torch.nn.ModuleList(PrimalConvBlocks)
        self.DualConvBlocks = torch.nn.ModuleList(DualConvBlocks)

        # Define forward and adjoint operators, handling necessary reshaping
        self.FwdOperator = lambda x: FwdOperator(x).squeeze(-3)  # Remove singleton dim
        self.AdjOperator = lambda x: AdjOperator(x.unsqueeze(-3).unsqueeze(-3)) # Add singleton dims

        self.relu = torch.nn.ReLU()

    def forward(self, y):
        """
        Performs a forward pass of the learned primal-dual algorithm.

        Args:
            y (torch.Tensor): The input sinogram data (acquisition data).

        Returns:
            torch.Tensor: The reconstructed image.
        """

        # Initialize dual variable (h) and primal variable (f) with zeros
        h = torch.zeros_like(y, device=y.device)
        dummy_f = self.AdjOperator(y)  # Create a dummy to get the correct shape
        f = torch.zeros_like(dummy_f, device=y.device)
        n_iter = len(self.PrimalConvBlocks)  # Number of iterations

        for i in range(n_iter):
            # Dual Update
            Op_f = self.FwdOperator(f)  # Apply forward operator to primal variable
            dual_input = torch.cat([Op_f, y], dim=1) # Concatenate with input sinogram
            h = self.DualConvBlocks[i](dual_input)   # Apply dual convolutional block
            h = self.relu(h)                       # Apply ReLU activation

            # Primal Update
            OpAdj_h = self.AdjOperator(h)   # Apply adjoint operator to dual variable
            primal_input = torch.cat([f, OpAdj_h], dim=1) # Concat current f and A^T h
            f = self.PrimalConvBlocks[i](primal_input)    # Apply primal convolutional block
            f = self.relu(f)                         # Apply ReLU activation

        return f



class ConvBlock(torch.nn.Module):
    """
    A basic convolutional block with three convolutional layers and PReLU activations.
    """
    def __init__(self, in_channels=2, out_channels=1, kernel_size=3, hidden_channels=10):
        """
        Initializes the ConvBlock.

        Args:
            in_channels (int): Number of input channels.
            out_channels (int): Number of output channels.
            kernel_size (int): Size of the convolutional kernel.
            hidden_channels (int): Number of channels in the hidden layers.
        """
        super(ConvBlock, self).__init__()
        self.conv1 = torch.nn.Conv2d(in_channels, hidden_channels, kernel_size, padding="same", bias=False)
        self.conv2 = torch.nn.Conv2d(hidden_channels, hidden_channels, kernel_size, padding="same", bias=False)
        self.conv3 = torch.nn.Conv2d(hidden_channels, out_channels, kernel_size, padding="same", bias=False)
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        """
        Performs a forward pass through the convolutional block.
        """
        x = self.relu(self.conv1(x))
        x = self.relu(self.conv2(x))
        x = self.relu(self.conv3(x))
        return x


class PETVarNet(torch.nn.Module):
    """
    A PyTorch module implementing a very basic variational network for PET image reconstruction.
    """
    def __init__(self, ConvBlocks, objfuncgrad):
        """
        Initializes the PETVarNet module.

        Args:
            ConvBlocks (list):  List of convolutional blocks.
            objfuncgrad (callable):  Gradient of the objective function.
        """
        super(PETVarNet, self).__init__()
        self.ConvBlocks = torch.nn.ModuleList(ConvBlocks)
        self.objfuncgrad = objfuncgrad  # Gradient of the objective function
        self.relu = torch.nn.ReLU()

    def forward(self, x):
        """
        Forward pass of the variational network.

        Args:
            x (torch.Tensor): Input image.

        Returns:
            torch.Tensor: Reconstructed image.
        """
        for ConvBlock in self.ConvBlocks:
            # Each iteration: Convolutional block output + gradient descent step
            x = self.relu(ConvBlock(x) + self.objfuncgrad(x))
        return x



class UseCases:
    """
    Class to demonstrate different use cases of the networks.
    """
    def __init__(self, modality, dim, test_name):
        """
        Initializes the UseCases class.

        Args:
            modality (str): Imaging modality ('PET').
            dim (str):  Dimensionality ('2D').
            test_name (str):  Name of the test to run ('all').
        """
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        if modality == 'PET':
            if dim == '2D':
                self.acq_data, self.image, self.acq_model = get_pet_2d()
                # Convert SIRF data to PyTorch tensors
                self.torch_acq_data = sirf_to_torch(self.acq_data, self.device)
                self.torch_image = sirf_to_torch(self.image, self.device).unsqueeze(0) # Add batch dimension
            else:
                raise ValueError("Only 2D PET is supported")
        else:  # Added this else for completeness
            raise ValueError("Only PET modality is supported")

        if test_name == 'all':
            if dim == '2D':
                try:
                    self.pet_varnet()
                except Exception as e:
                    print(f"Error encountered: {e}")
                    print("PET Varnet failed")

                self.learned_primal_dual()
                self.compare_gradient_descents()
        else:
            raise ValueError("Invalid test name. Only 'all' is supported")




    def learned_primal_dual(self):
        """
        Demonstrates the learned primal-dual network.
        """
        print("Learned Primal Dual")

        # Create primal and dual convolutional blocks
        PrimalConvBlocks = [ConvBlock() for i in range(3)]
        DualConvBlocks = [ConvBlock() for i in range(3)]

        # Set up forward and adjoint operators using SIRF's Operator class
        FwdOperator = Operator(self.acq_model, self.image.get_uniform_copy(1))
        AdjOperator = Operator(sirf.AdjointOperator(self.acq_model.get_linear_acquisition_model()), self.acq_data.get_uniform_copy(1))

        # Create the network and move it to the device
        net = PETLearnedPrimalDual(PrimalConvBlocks, DualConvBlocks, FwdOperator, AdjOperator)
        net.to(self.device)

        # Prepare input data and define the loss function
        torch_input = sirf_to_torch(self.acq_data, self.device)
        relu = torch.nn.ReLU()
        loss = torch.nn.MSELoss()  # Mean Squared Error loss

        # Set up the optimizer (Adam)
        optimizer = torch.optim.Adam(net.parameters(), lr=1e-4)

        # Training loop
        for i in range(50):
            optimizer.zero_grad()       # Clear gradients
            pred = relu(net(torch_input))  # Forward pass
            loss_val = loss(pred, self.torch_image)  # Calculate loss
            loss_val.backward()        # Backpropagate the loss
            optimizer.step()          # Update network parameters
            print("Iteration: ", i, "Loss: ", loss_val.item())

        # Display and save the reconstructed image
        out = relu(net(torch_input))  # Final forward pass after training
        plt.imshow(out.detach().cpu().numpy()[0, 0])  # Show the first image in the batch
        plt.savefig("learned_primal_dual.png")
        plt.close()


    def pet_varnet(self):
        """
        Demonstrates the PET variational network.
        """
        print("PET Varnet")

        # Create convolutional blocks
        ConvBlocks = [ConvBlock(in_channels=1) for i in range(5)]

        # Create the objective function and its gradient
        obj_fun = pet.make_Poisson_loglikelihood(self.acq_data)
        obj_fun.set_acquisition_model(self.acq_model)
        obj_fun.set_up(self.image.get_uniform_copy(1))
        objfuncgrad = ObjectiveFunctionGradient(obj_fun, self.image.get_uniform_copy(1))

        # Create the variational network
        net = PETVarNet(ConvBlocks, objfuncgrad)
        net.to(self.device)

        # Input:  Initialize with ones
        torch_input = torch.ones_like(self.torch_image)
        relu = torch.nn.ReLU()
        loss = torch.nn.MSELoss()

        optimizer = torch.optim.Adam(net.parameters(), lr=1e-4)
        for i in range(10):
            optimizer.zero_grad()
            pred = relu(net(torch_input))  # Apply the network
            loss_val = loss(pred, self.torch_image) # Compare with target
            loss_val.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss_val.item())

        out = relu(net(torch_input))
        plt.imshow(out.detach().cpu().numpy()[0, 0])
        plt.savefig("pet_varnet.png")
        plt.close()

    def compare_gradient_descents(self):
        """
        Compares two gradient descent methods: one using the acquisition model
        and the other using the objective function directly.
        """
        print("Comparing Gradient Descents")
        lr = 0.1
        n_iter = 5
        print(f"Learning Rate: {lr}, Number of Iterations: {n_iter}")

        print("Gradient Descent with Acquisition Model")
        out_acq_model = self.gradient_descent_with_acq_model(lr=lr, n_iter=n_iter)

        print("Gradient Descent with Objective Function")
        out_obj_func = self.gradient_descent_with_obj_func(lr=lr, n_iter=n_iter)

        plt.figure(figsize=(12, 4))


        plt.subplot(1, 3, 1)
        plt.imshow(self.image.get_uniform_copy(1).as_array()[0])
        plt.colorbar()
        plt.title('Initial Image')


        plt.subplot(1, 3, 2)
        plt.imshow(out_acq_model)
        plt.colorbar()
        plt.title('Acq model wrapped')


        plt.subplot(1, 3, 3)
        plt.imshow(out_obj_func)
        plt.colorbar()
        plt.title('Obj func wrapped')

        # Adjust layout and show plot
        plt.tight_layout()
        plt.show()
        plt.savefig("gd_comparison.png")
        plt.close()


        print("Sum of absolute differences between GD solns: ", np.abs(out_acq_model - out_obj_func).sum())


    def gradient_descent_with_acq_model(self, lr, n_iter):
        """
        Performs gradient descent using the acquisition model within the loss.
        """
        # Initialize image parameters (to be optimized)
        torch_image = sirf_to_torch(self.image.get_uniform_copy(1), self.device).unsqueeze(0) # Add batch dim.
        torch_image_params = torch.nn.Parameter(torch_image)

        torch_measurements = sirf_to_torch(self.acq_data, self.device)
        torch_acq_model = Operator(self.acq_model, self.image.get_uniform_copy(1))  # Wrap in Operator

        relu = torch.nn.ReLU()
        # Use Poisson NLL loss (appropriate for PET data)
        loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum')
        torch_obj_func = lambda x: loss(torch_acq_model(relu(x)), torch_measurements)  # Define loss function

        optimizer = torch.optim.Adam([torch_image_params], lr=lr)

        for i in range(n_iter):
            optimizer.zero_grad()
            loss_val = torch_obj_func(torch_image_params)  # Compute the loss
            loss_val.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss_val.item())

        return torch_image_params.data.detach().cpu().squeeze().numpy()  # Return optimized image

    def gradient_descent_with_obj_func(self, lr, n_iter):
        """
        Performs gradient descent directly on the SIRF objective function.
        """
        torch_image = sirf_to_torch(self.image.get_uniform_copy(1), self.device).unsqueeze(0)
        torch_image_params = torch.nn.Parameter(torch_image)  # Make it a parameter

        # Create the SIRF objective function
        obj_fun = pet.make_Poisson_loglikelihood(self.acq_data)
        obj_fun.set_acquisition_model(self.acq_model)
        obj_fun.set_up(self.image.get_uniform_copy(1))

        # Wrap the objective function for use with pytorch
        obj = ObjectiveFunction(obj_fun, self.image.get_uniform_copy(1))
        relu = torch.nn.ReLU()
        obj_relu = lambda x: -obj(relu(x)) # make negative for minimisation
        optimizer = torch.optim.Adam([torch_image_params], lr=lr)
        for i in range(n_iter):
            optimizer.zero_grad()
            loss = obj_relu(torch_image_params)  # Compute the loss
            loss.backward()
            optimizer.step()
            print("Iteration: ", i, "Loss: ", loss.item())
        return torch_image_params.data.detach().cpu().squeeze().numpy()  # Return optimised image



if __name__ == '__main__':
    torch.manual_seed(42)
    uses = UseCases('PET', '2D', 'all')
