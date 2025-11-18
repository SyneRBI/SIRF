# SIRF-PyTorch Wrapper
This wrapper provides a bridge between the [SIRF](https://github.com/SyneRBI/SIRF) (Synergistic Image Reconstruction Framework) library and [PyTorch](https://github.com/pytorch/pytorch), enabling the use of SIRF's image reconstruction operators and objective functions within PyTorch's automatic differentiation (autodiff) framework.

## Usage and Use Cases

The `sirf.torch.Operator`, `sirf.torch.ObjectiveFunction`, and `sirf.torch.ObjectiveFunctionGradient` classes are designed to be used as standard PyTorch `nn.Module`s.  You would initialise them with the appropriate SIRF objects and then use them in your forward pass like any other PyTorch layer.

`tests/use_cases.py` demonstrates `sirf.torch` integration in PyTorch with minimal 2D PET examples:

*   **Learned Primal-Dual:** Implements a learned primal-dual network for PET image reconstruction, showcasing the use of `sirf.torch.Operator` for handling the forward and adjoint projection operations.
*   **PET Variational Network (PETVarNet):**  Demonstrates a variational network approach, combining convolutional blocks with gradient information from a SIRF objective function using `sirf.torch.ObjectiveFunctionGradient`.
*   **ADAM Gradient Descent Comparison:**  Compares two gradient descent implementations: one leveraging the `sirf.torch.Operator` for the acquisition model within the loss calculation, and another directly utilising the `sirf.torch.ObjectiveFunction` for a more traditional optimisation approach.  This highlights the flexibility of the wrapper in different optimisation strategies.

## Dimensions used by the wrapper

The wrappers prioritise SIRF's data formats, meaning that the torch arrays must have the shape:
* [batch, [channel,] *SIRF.DataContainer.shape], where the channel dimension is optional.

This requires the **user** to ensure the dimensionality to match between layers.

### Example dimension manipulation
For example, a sinogram in SIRF has shape [tof bins, sinograms, views, tang pos]. For a single non-tof sinogram this is [1, 1, views, tang pos]. The expected torch tensor shape for this wrapper is [batch, [channel,] 1, 1, views, tang pos]. On the otherhand a 2D convolution requires [batch, [channel,] height, width].

```python
conv_1 = torch.nn.Conv2D()
conv_2 = torch.nn.Conv2D()
adjoint_operator = sirf.SIRF_torch.Operator(sirf.AdjointOperator(acquisition_model))
y # sinogram of dimension [batch, [channel,] views, tang pos]

y_filtered = conv_1(y) # filtered sinogram of dimension [batch, [channel,] views, tang pos]
y_filtered = y_filtered.unsqueeze(-3).unsqueeze(-3) # filtered sinogram of dimension [batch, [channel,] 1, 1, views, tang pos]
x_bp = adjoint_operator(y_filtered) # back-projected image of dimension [batch, [channel,] 1, height, width]
x_bp = x_bp.squeeze(-3)  # back-projected image of dimension [batch, [channel,] height, width]
x_bp_filtered = conv_2(x_bp)  # filtered back-projected image of dimension [batch, [channel,] height, width]
```

## Forward and backward clarification

The use of the terms forward and backward have different meaning given the context:
* Automatic differentiation: Forward (tangent) mode autodiff computes the Jacobian-Vector-Product (JVP). This propagates derivatives forward along with the function evaluation. Backward (or reverse/adjoint) mode autodiff is the Vector-Jacobian-Product (VJP) that propagates derivative information in the reverse direction of the function's evaluation. 
* Backward autodiff cont.: Forward pass evaluates the function saving intermediate values. Backward pass uses the chain rule and intermediate values computing the derivatives in the reverse direction with the VJP.
* `torch.autograd.Function`: the `forward` method (forward pass) is the function evaluation. The `backward` method (backward pass) computes the VJP. More specifically, the `backward(*grad_output)` method multiplies the `grad_output` which represents the gradient(s) of a subsequent function/operator (evaluated at the output of `forward`), via chain-rule, by the adjoint of the Jacobian of the `forward` method. 

This SIRF-PyTorch wrapper is **only** for reverse-mode automatic differentiation via subclassing `torch.autograd.Function`.

## Wrapper Design

The wrapper provides three main classes:

1.  `sirf.torch.Operator`: Wraps a SIRF `Operator` (e.g., a projection operator). Applies the operator forward pass, and applies the adjoint of the Jacobian in backward pass.
2.  `sirf.torch.ObjectiveFunction`: Wraps a SIRF `ObjectiveFunction` for computing its value in the forward pass, and multiplying with the objective function gradient in the backward pass.
3.  `sirf.torch.ObjectiveFunctionGradient`: Wraps a SIRF `ObjectiveFunction` that computes the objective function gradients in the forward pass and the Hessian-vector product in the backward pass. In the backward the Hessian is evaluated at the point which the objective function's gradient was evaluated.

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
    1.  Converts the input PyTorch tensor (representing an image for instance) to a SIRF object.
    2.  Calls the SIRF `ObjectiveFunction.__call__()` method.
    3.  Returns the of the objective function value as a PyTorch tensor.
    4.  Saves relevant information to the `ctx` if gradients are needed.

*   **Backward Pass (VJP):**
    1.  Receives the upstream gradient (`grad_output`), in this case it is always a scalar.
    2.  Gets the gradient of the objective function using `sirf_obj_func.gradient()`, which computed at the input and multiplied by the upstream gradient.
    3.  Converts the SIRF gradient to a PyTorch tensor.
    4.  Returns the gradient multiplied by `grad_output`.


### `_ObjectiveFunctionGradient` (Forward and Backward Passes)

*   **Forward Pass:**
    1.  Converts the input PyTorch tensor to a SIRF object.
    2.  Computes the *gradient* of the SIRF objective function using `sirf_obj_func.gradient()`, which is computed on the input.
    3.  Returns the gradient as a PyTorch tensor.

*   **Backward Pass (VJP):**
    1.  Receives the upstream gradient (`grad_output`), which now represents a *vector* (not a scalar) of the same shape as the output of `forward`.
    2.  Converts `grad_output` to a SIRF object.
    3.  Multiples the Hessian evaluated at the input of `forward` with the "upstream gradient" using `sirf_obj_func.multiply_with_Hessian()`.
    4.  Returns the Hessian multiplied with a vector as a tensor.

# TODO

* Extend to subsets in the wrapper
* Extend objective functions that vary between batch items
