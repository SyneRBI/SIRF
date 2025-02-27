try:
    import torch
    import sirf
    
    from sirf.STIR import *
except ModuleNotFoundError:
    raise ModuleNotFoundError('Failed to import torch. Please install PyTorch first.')

# based on 
# https://github.com/educating-dip/pet_deep_image_prior/blob/main/src/deep_image_prior/torch_wrapper.py


def sirf_to_torch(
        sirf_src: object,
        device: torch.device,
        requires_grad: bool = False
        ) -> torch.Tensor:

    # use torch.tensor to infer data type
    return torch.tensor(sirf_src.as_array(), requires_grad=True).to(device)

def torch_to_sirf_(
        torch_src: torch.Tensor,
        sirf_dest: object,
        ) -> object:

    # This is an in-place operation - CAREFUL
    # Only really to be used within torch.autograd.Function
    return sirf_dest.fill(torch_src.detach().cpu().numpy())


class _ObjectiveFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image: torch.Tensor,
            sirf_image_template: object,
            sirf_obj_func: object
            ):

        device = torch_image.device
        sirf_image_template = torch_to_sirf_(torch_image, sirf_image_template)
        # Negative for Gradient Descent as per torch convention
        # HERE WE ARE COPYING THE VALUE
        value = torch.tensor(-sirf_obj_func.get_value(sirf_image_template)).to(device)
        if torch_image.requires_grad:
            ctx.device = device
            ctx.sirf_image_template = sirf_image_template
            ctx.sirf_obj_func = sirf_obj_func
            # ensure value is a tensor with requires_grad=True
            return value.requires_grad_(True)
        else:
            return value

    @staticmethod
    def backward(ctx,
            grad_output: torch.Tensor
            ):
            
        sirf_obj = ctx.sirf_obj_func
        sirf_image_template = ctx.sirf_image_template
        device = ctx.device
        # Negative for Gradient Descent as per torch convention
        sirf_grad = -sirf_obj.get_gradient(sirf_image_template)
        grad = sirf_to_torch(sirf_grad, device, requires_grad=True)
        return grad_output*grad, None, None, None


class _ObjectiveGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, obj_fnc):
        # We consider the objective function f(x): R^n -> R (maps n-dimensional input to a scalar).
        # This forward pass computes the gradient of f(x) with respect to x.
        # g(x) = ∇f(x) = J^T  (where J is the Jacobian of f(x)).  g(x) is a row vector.
        # The VJP is [1] * J = ∇f(x), where [1] is a scalar (or a 1x1 tensor).
        # We use [1] because we're computing the full gradient of the scalar output f(x).
        # The result of the VJP is the gradient itself.  The "upstream gradient"
        # is implicitly [1] because there are no prior operations.
        # Here we pass this gradient.
        raise NotImplementedError

    @staticmethod
    def backward(ctx, grad_output):
        # Compute the Hessian-vector product (HVP).
        # Explanation:
        # v = grad_output = ∂L/∂g(x), where L is the loss function and g(x) is the
        # *conceptual* output of this Function (the gradient of f(x)).  v is a vector.
        # The Hessian is ∇²f(x) = H|x, evaluated at the input point x.
        # The HVP is v^T * H|x.
        # This HVP computes ∂L/∂x, the gradient of the loss with respect to the *input* (x)
        # of the forward pass, by applying the chain rule: ∂L/∂x = (∂L/∂g(x)) * (∂g(x)/∂x) = v^T * H.
        raise NotImplementedError


class _ObjectiveGradient(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, obj_fnc):
        # We consider the objective function f(x): R^n -> R (maps n-dimensional input to a scalar).
        # This forward pass represent* the computation of the Vector-Jacobian Product (VJP).
        # g(x) = ∇f(x) = J^T  (where J is the Jacobian of f(x) at point x).  g(x) is a row vector.
        # The VJP, in general, is v^T * J.  Within this forward pass, we are
        # computing the full gradient, so v = [1] (a scalar or 1x1 tensor).
        # The result of this VJP *is* the gradient itself: [1] * J = ∇f(x).
        raise NotImplementedError

    @staticmethod
    def backward(ctx, grad_output):
        # This backward pass represents the computation of the Hessian-Vector Product (HVP).
        # Explanation:
        # v = grad_output = ∂L/∂g(x), where L is the loss function and g(x) is the
        # *conceptual* output of the forward pass (the gradient of f(x)). v is a vector.
        # The Hessian, H, is ∇²f(x) = H|x, evaluated at the input point x.
        # The HVP is v^T * H|x.
        # This HVP computes ∂L/∂x, the gradient of the loss with respect to the *input* (x)
        # of the forward pass, by applying the chain rule: ∂L/∂x = (∂L/∂g(x)) * (∂g(x)/∂x) = v^T * H.
        raise NotImplementedError


class _Operator(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_src: torch.Tensor,
            sirf_src_template: object,
            sirf_operator: object
            ):

        device = torch_src.device
        sirf_src_template = torch_to_sirf_(torch_src, sirf_src_template)
        sirf_dest = sirf_operator.forward(sirf_src_template)
        if torch_src.requires_grad:
            ctx.device = device
            ctx.sirf_dest = sirf_dest
            ctx.sirf_operator = sirf_operator
            return sirf_to_torch(sirf_dest, device, requires_grad=True)
        else:
            return sirf_to_torch(sirf_dest, device)

    @staticmethod
    def backward(ctx,
            grad_output: torch.Tensor
            ):

        sirf_src = ctx.sirf_operator.backward(torch_to_sirf_(grad_output, ctx.sirf_dest))
        grad = sirf_to_torch(sirf_src, ctx.device, requires_grad=True)
        return grad, None, None, None

def check_shapes(torch_shape, sirf_shape):
    if torch_shape != sirf_shape:
        raise ValueError(f"Invalid shape. Expected sirf shape {sirf_shape} but \
            got torch shape {torch_shape}")

class SIRFOperator(torch.nn.Module):
    def __init__(self,
            operator, 
            sirf_src_template
            ):
        super(SIRFOperator, self).__init__()
        # get the shape of src
        self.operator = operator
        self.sirf_src_shape = sirf_src_template.as_array().shape
        self.sirf_src_template = sirf_src_template
        self.channel = False
        
    def forward(self, torch_src):
        # PyTorch src is size [batch, channel, *src_shape] or
        #  [batch, *src_shape]

        # check the src shape
        torch_src_shape = torch_src.shape
        if len(torch_src_shape) == len(self.sirf_src_shape):
            raise ValueError(f"Invalid shape of src. Expected batch dim.")
        elif len(torch_src_shape) == len(self.sirf_src_shape) + 1:
            check_shapes(torch_src_shape[1:], self.sirf_src_shape)
            if self.sirf_src_shape == torch_src_shape[1:]:
                torch_src = torch_src.unsqueeze(1) # add channel dimension
                self.channel = True
        elif len(torch_src_shape) == self.sirf_src_shape + 2:
            check_shapes(torch_src_shape[2:], self.sirf_src_shape)

        n_batch = torch_src.shape[0]
        n_channel = torch_src.shape[1]

        # This looks horrible, but PyTorch will be able to trace.
        batch_dest = []
        for batch in range(n_batch):
            channel_dest = []
            for channel in range(n_channel):
                channel_dest.append(_Operator.apply(torch_src[batch, channel], 
                    self.sirf_src_template, self.operator))
            batch_dest.append(torch.stack(channel_dest, dim=0))
        if self.channel:
            return torch.stack(batch_dest, dim=0).squeeze(1)
        # [batch, channel, *sirf_dest.shape]
        return torch.stack(batch_dest, dim=0)


class SIRFObjectiveFunction(torch.nn.Module):
    def __init__(self,
            sirf_obj_func: object,
            sirf_image_template: object
            ):
        super(SIRFObjectiveFunction, self).__init__()
        self.sirf_obj_func = sirf_obj_func
        self.sirf_image_template = sirf_image_template
        self.sirf_image_shape = sirf_image_template.shape

    def forward(self, torch_image):
        # PyTorch src is size [batch, channel, *sirf_image_shape] or
        #  [batch, *sirf_image_shape]

        torch_image_shape = torch_image.shape
        if len(torch_image_shape) == len(self.sirf_image_shape):
            raise ValueError(f"Invalid shape of src. Expected batch dim.")
        elif len(torch_image_shape) == len(self.sirf_image_shape) + 1:
            check_shapes(torch_image_shape[1:], self.sirf_image_shape)
            if self.sirf_image_shape == torch_image_shape[1:]:
                torch_image = torch_image.unsqueeze(1) # add channel dimension
                self.channel = True
        elif len(torch_image_shape) == len(self.sirf_image_shape) + 2:
            check_shapes(torch_image_shape[2:], self.sirf_image_shape)

        n_batch = torch_image.shape[0]
        n_channel = torch_image.shape[1]

        # This looks horrible, but PyTorch should be able to trace.
        batch_values = []
        for batch in range(n_batch):
            channel_values = []
            for channel in range(n_channel):
                channel_values.append(_ObjectiveFunction.apply(torch_image[batch, channel], self.sirf_image_template, self.sirf_obj_func))
            batch_values.append(torch.stack(channel_values, dim=0))
        # [batch, channel, *forward_projected.shape]
        return torch.stack(batch_values, dim=0) 

