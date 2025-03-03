# SIRF-PyTorch Wrapper
This wrapper provides a bridge between the [SIRF](https://github.com/SyneRBI/SIRF) (Synergistic Image Reconstruction Framework) library and PyTorch, enabling the use of SIRF's image reconstruction operators and objective functions within PyTorch's automatic differentiation (autodiff) framework. 

## Core Concepts: Autodiff, Reverse Mode, JVP, VJP, and HVP

This wrapper leverages PyTorch's `torch.autograd.Function` to achieve integration with SIRF.  Understanding the following concepts is crucial for the design:

*   **Automatic Differentiation (Autodiff):** A technique for numerically evaluating the derivative of a function specified by a computer program.  PyTorch uses a "reverse mode" autodiff system.

*   **Reverse Mode Autodiff:**  In reverse mode, the gradients are computed by traversing the computational graph backwards, from the output (typically a scalar loss function) to the inputs.  This is highly efficient for functions with many inputs and a single output, which is common in machine learning.

*   **Jacobian-Vector Product (JVP):**  Let `y = f(x)` be a vector-valued function, where `x` is an input vector and `y` is an output vector. The Jacobian matrix `J` contains all the first-order partial derivatives of `f`.  The JVP is the product `Jv`, where `v` is a vector.  It gives the directional derivative of `f` at `x` in the direction of `v`.  *Forward mode* autodiff computes JVPs.

*   **Vector-Jacobian Product (VJP):** With `y = f(x)` as above, the VJP is the product `w^T J`, where `w` is a vector (often called the "upstream/puput gradient").  The VJP represents the gradient of a scalar-valued function that depends on `y` (e.g., a loss function), with respect to the input `x`.  *Reverse mode* autodiff computes VJPs. 

*   **Hessian-Vector Product (HVP):** The Hessian matrix `H` contains all the second-order partial derivatives of a scalar-valued function `f(x)`. The HVP is the product `Hv`, where `v` is a vector.  It gives the second-order directional derivative of `f` at `x` in the direction of `v`. 

## Wrapper Design

The wrapper provides three main classes:

1.  `SIRFTorchOperator`: Wraps a SIRF `Operator` (e.g., a projection operator). Computes the JVP forward and VJP in reverse.
2.  `SIRFTorchObjectiveFunction`: Wraps a SIRF `ObjectiveFunction` for calculating its value forward, and gradient in reverse mode.
3.  `SIRFTorchObjectiveFunctionGradient`: Wraps a SIRF `ObjectiveFunction` for calculating its gradient forward and computing the Hessian-vector product in reverse mode.

These classes use custom `torch.autograd.Function` implementations (`_Operator`, `_ObjectiveFunction`, and `_ObjectiveFunctionGradient`) to define the forward and backward passes, handling the conversions between PyTorch tensors and SIRF objects.

### `_Operator` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor to a SIRF object.
    2.  Applies the SIRF `Operator.forward()` method.  Let's represent the SIRF Operator as `A`.  So, this step computes `y = Ax`, where `x` is the input (SIRF object) and `y` is the output (SIRF object).
    3.  Converts the result back to a PyTorch tensor.
    4.  If the input tensor requires gradients, it saves relevant information (the output SIRF object and the operator) in the context (`ctx`) for use in the backward pass.

*   **Backward Pass (Adjoint Operation - VJP):**
    1.  Receives the "upstream gradient" (`grad_output`) – the gradient of the loss `L` with respect to the *output* of the forward pass: `∂L/∂y`.
    2.  Converts `grad_output` to a SIRF object.
    3.  Applies the SIRF `Operator.backward()` method (which computes the adjoint operation). The adjoint operation, `A^T`, effectively computes the VJP.  If the forward pass computes `y = Ax`, the backward pass computes `x̄ = A^T ȳ`, where `A^T` is the adjoint of `A`, `ȳ` is `grad_output = ∂L/∂y`, and `x̄` is the gradient with respect to the input `x`. So, `x̄ = ∂L/∂x = (∂L/∂y)(∂y/∂x) =  ȳ^T J = (A^T)ȳ`. Here, `J` is the Jacobian of the forward operation `y=Ax`, and it is implicitly represented by the operator `A` itself.
    4.  Converts the resulting SIRF object back to a PyTorch tensor and returns it. This is the gradient with respect to the *input* of the forward pass (`∂L/∂x`).

    **Connection to VJP:** The backward pass directly computes the VJP.  The upstream gradient `grad_output` (which is `∂L/∂y`) acts as the vector `w^T` in the VJP formula `w^T J`. The SIRF `Operator.backward()` method implicitly performs the multiplication by the Jacobian transpose (`J^T`, which is equivalent to the adjoint `A^T`).

### `_ObjectiveFunction` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor (representing an image) to a SIRF object.
    2.  Calls the SIRF `ObjectiveFunction.get_value()` method. Let `f(x)` represent the SIRF objective function. This step computes `f(x)`.
    3.  Returns the *negative* of the objective function value as a PyTorch tensor: `-f(x)`. We negate the value because PyTorch optimizers perform *minimization*, and we typically want to maximize an objective function (or minimize its negative) in image reconstruction.
    4. Saves relevant information to the `ctx` if gradients are needed.

*   **Backward Pass (VJP):**
    1.  Receives the upstream gradient (`grad_output`), which is typically a tensor containing the scalar `1` (since the output of the forward pass is a scalar, and `∂L/∂(-f(x)) = 1` if `L = -f(x)`).
    2.  Gets the gradient of the *negative* objective function using `sirf_obj_func.get_gradient()`, which computes `-∇f(x)`.
    3.  Converts the SIRF gradient to a PyTorch tensor.
    4.  Returns the gradient multiplied by `grad_output`.

    **More on the VJP:** Let `f(x)` be the SIRF objective function. The forward pass computes `-f(x)`.  The backward pass needs to compute `∂L/∂x`, where `L` is the overall loss function (which, in the simplest case, *is* the output of the forward pass, `-f(x)`). Using the chain rule:
    `∂L/∂x = (∂L/∂(-f(x))) * (∂(-f(x))/∂x) = grad_output * (-∇f(x))`.  Since the forward pass output is a scalar, the Jacobian is simply the gradient `∇f(x)`.  The upstream gradient `grad_output` represents `∂L/∂(-f(x))`, and `sirf_obj_func.get_gradient()` provides `-∇f(x)`.  The multiplication `grad_output * (-∇f(x))` correctly computes the VJP.

### `_ObjectiveFunctionGradient` (Forward and Backward Passes)

*   **Forward Pass (JVP with v=1):**
    1.  Converts the input PyTorch tensor to a SIRF object.
    2.  Computes the *gradient* of the *negative* SIRF objective function using `sirf_obj_func.get_gradient()`, which returns `-∇f(x)`.  This is equivalent to computing a JVP where the input vector `v` is implicitly `1` (because we are taking the full gradient of the *scalar* output of the objective function). The result is a vector, the gradient.
    3.  Returns the (negative) gradient as a PyTorch tensor.

*   **Backward Pass (Hessian-Vector Product - HVP):**
    1.  Receives the upstream gradient (`grad_output`), which now represents a *vector* (not a scalar).  This vector corresponds to `∂L/∂g(x)`, where `g(x) = -∇f(x)` is the *output* of the forward pass (the negative gradient of the objective function).
    2.  Converts `grad_output` to a SIRF object.
    3.  Computes the Hessian-vector product (HVP) using `sirf_obj_func.multiply_with_Hessian()`.  This method takes two arguments: the current estimate (`sirf_image`, where the Hessian is evaluated) and the input vector (`sirf_grad`, which is `grad_output` converted to a SIRF object).
    4. Returns the (negative) HVP.

    **More on the HVP:**  Let `f(x)` be the SIRF objective function, and let `g(x) = -∇f(x)` be the output of the forward pass. We want to compute `∂L/∂x`, where `L` is the overall loss. The backward pass receives `v = grad_output = ∂L/∂g(x)`.  The Hessian of `f(x)` is `H = ∇²f(x)`. We compute the HVP: `H|x * v`, where `H|x` denotes the Hessian evaluated at point x. Applying the chain rule:

    `∂L/∂x = (∂L/∂g(x)) * (∂g(x)/∂x) = v^T * (∂(-∇f(x))/∂x) = - v^T * ∇²f(x) = -v^T * H|x`.

    The result of `sirf_obj_func.multiply_with_Hessian` is negated in the code, consistent with the negative gradient used throughout. This ensures we're performing gradient *descent*.

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
loss = torch.nn.PoissonNLLLoss(log_input=False, full=False) # Example loss
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
* The negative introduced for PET is not neccessary for MR/CIL
