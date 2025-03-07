# SIRF-PyTorch Wrapper
This wrapper provides a bridge between the [SIRF](https://github.com/SyneRBI/SIRF) (Synergistic Image Reconstruction Framework) library and PyTorch, enabling the use of SIRF's image reconstruction operators and objective functions within PyTorch's automatic differentiation (autodiff) framework. 

## Forward and backward clarification

The use of the terms forward and backward have different meaning given the context:
* Automatic differentiation: Forward (tangent) mode autodiff computes the Jacobian-Vector-Product (JVP). This propagates derivatives forward along with the function evaluation. Backward (or reverse/adjoint) mode autodiff is the Vector-Jacobian-Product (VJP) that propagates derivative information in the reverse direction of the function's evaluation. 
* Reverse-mode autodiff: Forward pass evaluates the function saving intermediate values. Backward pass uses the chain rule and intermediate values computing the derivatives in the reverse direction with the VJP.
* `torch.autograd.Function`: the `forward` method (forward pass) is the function evaluation. The `backward` method (backward pass) computes the VJP. More specifically, the `backward(*grad_output)` method multiplies the `grad_output` which represents the gradient(s) of a subsequent scalar-valued objective function with respect to the output of `forward`, via chain-rule, by the adjoint of the Jacobian of the `forward` method. 

This SIRF-PyTorch wrapper is **only** for reverse-mode automatic differentiation via subclassing `torch.autograd.Function`.

## Wrapper Design

The wrapper provides three main classes:

1.  `SIRFTorchOperator`: Wraps a SIRF `Operator` (e.g., a projection operator). Applies the operator forward pass, and applies the adjoint in backward pass.
2.  `SIRFTorchObjectiveFunction`: Wraps a SIRF `ObjectiveFunction` for computes the value in the forward pass, and objective function gradient in backward pass.
3.  `SIRFTorchObjectiveFunctionGradient`: Wraps a SIRF `ObjectiveFunction` that computes the objective function gradients in the forward pass and the Hessian-vector product in the backward pass. In the backward the Hessian is evaluated at the point which the objective functions gradient was evaluated.

These classes use custom `torch.autograd.Function` implementations (`_Operator`, `_ObjectiveFunction`, and `_ObjectiveFunctionGradient`) to define the forward and backward passes, handling the conversions between PyTorch tensors and SIRF objects.

### `_Operator` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor to a SIRF object.
    2.  Applies the SIRF `Operator.forward()` method.
    3.  Converts the result back to a PyTorch tensor.
    4.  If the input tensor requires gradients, it saves relevant information (the output SIRF object and the operator) in the context (`ctx`) for use in the backward pass.

*   **Backward Pass (VJP):**
    1.  Receives the "upstream gradient" (`grad_output`).
    2.  Converts `grad_output` to a SIRF object.
    3.  Applies the SIRF `Operator.backward()` method. This will apply the **Jacobian adjoint** of the operator to upstream gradient (the vector).
    4.  Converts the resulting SIRF object back to a PyTorch tensor and returns it.

### `_ObjectiveFunction` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor (representing an image) to a SIRF object.
    2.  Calls the SIRF `ObjectiveFunction.get_value()` method.
    3.  Returns the *negative* of the objective function value as a PyTorch tensor. We negate the value because PyTorch optimizers perform *minimization*, and we typically want to maximize an objective function (or minimize its negative) in image reconstruction.
    4. Saves relevant information to the `ctx` if gradients are needed.

*   **Backward Pass (VJP):**
    1.  Receives the upstream gradient (`grad_output`), in this case it is always a scalar.
    2.  Gets the gradient of the *negative* objective function using `sirf_obj_func.get_gradient()`, which computed at the input and multiplied by the upstream gradient.
    3.  Converts the SIRF gradient to a PyTorch tensor.
    4.  Returns the gradient multiplied by `grad_output`.


### `_ObjectiveFunctionGradient` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor to a SIRF object.
    2.  Computes the *gradient* of the *negative* SIRF objective function using `sirf_obj_func.get_gradient()`, which is computed on the input.
    3.  Returns the (negative) gradient as a PyTorch tensor.

*   **Backward Pass (VJP):**
    1.  Receives the upstream gradient (`grad_output`), which now represents a *vector* (not a scalar).
    2.  Converts `grad_output` to a SIRF object.
    3.  Multiples the Hessian evaluated as the input with the upstream gradient using `sirf_obj_func.multiply_with_Hessian()`.
    4. Returns the Hessian multiplied with a vector.

## Usage and Use Cases

The `SIRFTorchOperator`, `SIRFTorchObjectiveFunction`, and `SIRFTorchObjectiveFunctionGradient` classes are designed to be used as standard PyTorch `nn.Module`s.  You would initialise them with the appropriate SIRF objects and then use them in your forward pass like any other PyTorch layer.

The following use cases illustrate how to leverage this wrapper for different image reconstruction scenarios.

### 1. Gradient Descent with Acquisition Model

This example demonstrates how to perform gradient descent using the wrapped SIRF acquisition model directly within the objective function. This is a basic approach, useful for simple reconstructions.

```python
import torch
from sirf.SIRF_torch import sirf_to_torch, SIRFTorchOperator
# Example with a SIRF Projector and a PyTorch loss function:

# Assuming:
# - acq_data, image, acq_model are pre-defined SIRF objects.
# - device is a torch.device (e.g., 'cuda' or 'cpu')

torch_image = sirf_to_torch(image.get_uniform_copy(1), device).unsqueeze(0)  # Initial image
torch_image_params = torch.nn.Parameter(torch_image) # Make it a trainable parameter
torch_measurements = sirf_to_torch(acq_data, device)
torch_acq_model = SIRFTorchOperator(acq_model, image.get_uniform_copy(1))
relu = torch.nn.ReLU() # Apply non-negativity
loss = torch.nn.PoissonNLLLoss(log_input=False, full=False, size_average=None, eps=1e-08, reduce=None, reduction='sum') # Example loss
torch_obj_func = lambda x: loss(torch_acq_model(relu(x)), torch_measurements)

optimizer = torch.optim.Adam([torch_image_params], lr=0.1)

for i in range(100):  # Number of iterations
    optimizer.zero_grad()
    loss_val = torch_obj_func(torch_image_params)
    loss_val.backward()
    optimizer.step()
    print(f"Iteration {i}, Loss: {loss_val.item()}")

# Reconstructed image is now in torch_image_params.data
```

# TODO

* Extend to subsets in the wrapper
* Extend objective functions that vary between batch items
* The negative introduced for PET is not neccessary for MR/CIL - but what is a CIL/MR objective? Could be better to focus on the acquisition models.
* Should the use cases just be the exercises?
