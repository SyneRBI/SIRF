
from sirf.torch import Operator, ObjectiveFunction, ObjectiveFunctionGradient, sirf_to_torch
import sirf.STIR as pet
pet.set_verbosity(1)
pet.AcquisitionData.set_storage_scheme("memory")
import sirf
msg = sirf.STIR.MessageRedirector(info=None, warn=None, errr=None)
from sirf.Utilities import examples_data_path
import matplotlib.pyplot as plt
import torch

pet_data_path = examples_data_path('PET')
pet_2d_raw_data_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/template_sinogram.hs')
pet_2d_acq_data = pet.AcquisitionData(pet_2d_raw_data_file)
pet_2d_init_image_file = pet.existing_filepath(pet_data_path, 'thorax_single_slice/emission.hv')
pet_2d_image_data = pet.ImageData(pet_2d_init_image_file)
pet_2d_acq_model = pet.AcquisitionModelUsingParallelproj()
data_processor = pet.TruncateToCylinderProcessor()
pet_2d_acq_model.set_image_data_processor(data_processor)
pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)
sens_img = pet_2d_acq_model.backward(pet_2d_acq_data.get_uniform_copy(1.0))
inv_sens_img = sens_img.power(-1)
data_processor.apply(inv_sens_img)
pet_2d_acq_data = pet_2d_acq_model.forward(pet_2d_image_data) + 10.5
pet_2d_acq_model.set_background_term(pet_2d_acq_data.get_uniform_copy(10.5))
pet_2d_acq_model.set_up(pet_2d_acq_data, pet_2d_image_data)

acq_data = pet_2d_acq_data
acq_model = pet_2d_acq_model
image_data = pet_2d_image_data
obj_fun = pet.make_Poisson_loglikelihood(acq_data)
obj_fun.set_acquisition_model(acq_model)
obj_fun.set_up(image_data)
objfuncgrad = ObjectiveFunctionGradient(obj_fun, image_data)

dev = "cuda" if torch.cuda.is_available() else "cpu"

cnn = torch.nn.Sequential(
    torch.nn.Conv2d(1, 5, 5, padding="same", bias=False),
    torch.nn.Conv2d(5, 5, 5, padding="same", bias=False),
    torch.nn.PReLU(device=dev),
    torch.nn.Conv2d(5, 5, 5, padding="same", bias=False),
    torch.nn.Conv2d(5, 5, 5, padding="same", bias=False),
    torch.nn.PReLU(device=dev),
    torch.nn.Conv2d(5, 1, 1, padding="same", bias=False),
).to(dev)


class UnrolledOSEMVarNet(torch.nn.Module):
    def __init__(
        self,
        objective_function_gradient: sirf.torch.ObjectiveFunctionGradient,
        inv_sens_img: torch.Tensor,
        convnet: torch.nn.Module,
        device: str,
    ) -> None:
        """Unrolled OSEM Variational Network with 2 blocks

        Parameters
        ----------
        objective_function : sirf.STIR objetive function
            (listmode) Poisson logL objective function
            that we use for the OSEM updates
        sirf_template_image : sirf.STIR.ImageData
            used for the conversion between torch tensors and sirf images
        convnet : torch.nn.Module
            a (convolutional) neural network that maps a minibatch tensor 
            of shape [1,1,spatial_dimensions] onto a minibatch tensor of the same shape
        device : str
            device used for the calculations
        """
        super().__init__()

        # OSEM update layer using the 1st subset of the listmode data
        self.objective_function_gradient = objective_function_gradient
        self._inv_sens_img = inv_sens_img

        self._convnet = convnet
        self._relu = torch.nn.ReLU()

        # trainable parameters for the fusion of the OSEM update and the CNN output in the two blocks
        # we start with a weight of 10 for the fusion
        # a good starting value depends on the scale of the input image
        self._fusion_weight0 = torch.nn.Parameter(
            10 * torch.ones(1, device=device, dtype=torch.float32)
        )
        self._fusion_weight1 = torch.nn.Parameter(
            10 * torch.ones(1, device=device, dtype=torch.float32)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x1 = self._relu(
            self._fusion_weight0 * self._convnet(x) + self.objective_function_gradient(x)*x*self._inv_sens_img
        )
        x2 = self._relu(
            self._fusion_weight1 * self._convnet(x1) + self.objective_function_gradient(x1)*x1*self._inv_sens_img
        )

        return x2

inv_sens_img = torch.tensor(inv_sens_img.as_array()).unsqueeze(0).to(dev)
varnet = UnrolledOSEMVarNet(objfuncgrad, inv_sens_img, cnn, dev)
varnet.to(dev)

# use seed for reproducibility
torch.manual_seed(42)

torch_image = torch.tensor(image_data.as_array()).unsqueeze(0).to(dev)
# add gaussian noise to the image
torch_input = torch_image + 10 * torch.randn_like(torch_image)
# clip the input to be non-negative
torch_input = torch.clamp(torch_input, 0).detach()
# plot non noisy image and noisy image
plt.imshow(torch_image.detach().cpu().numpy()[0,0])
plt.colorbar()
plt.savefig("pet_image.png")
plt.close()
plt.imshow(torch_input.detach().cpu().numpy()[0,0])
plt.colorbar()
plt.savefig("pet_noisy_image.png")
plt.close()

# set up the optimizer
optimizer = torch.optim.Adam(varnet.parameters(), lr=1e-4)
relu = torch.nn.ReLU()
loss = torch.nn.MSELoss()
for i in range(100):
    optimizer.zero_grad()
    loss_val = loss(relu(varnet(torch_input)),torch_image)
    loss_val.backward()
    optimizer.step()
    print("Iteration: ", i, "Loss: ", loss_val.item())
out = relu(varnet(torch_input))
plt.imshow(out.detach().cpu().numpy()[0,0])
plt.savefig("pet_varnet.png")
