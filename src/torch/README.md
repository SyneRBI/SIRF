# Notes for Design Discussion


## File structure

Currently the source is found here:
```md
SIRF
├── src
│   ├── torch
│       ├── CMakeLists.txt
│       ├── README.md
│       └── SIRF_torch.py
```
Currently installed in:
```md
INSTALL
├── python
│   ├── sirf
│       └── SIRF_torch.py
```

The classes are accessed as `sirf.SIRF_torch.TheClass`.

## What is within `SIRF_torch`?

Briefly going from the most abstracted:

---

So we have two functions for transferring from sirf to torch, and back. The reason for abstracting is because it is used often, as well as being somewhat prudent if this changes in the future. For example if we don't need to copy anymore.

```python
def sirf_to_torch(sirf_src, torch_dest, requires_grad=False):
    if requires_grad:
        return torch.tensor(sirf_src.as_array(), requires_grad=True).to(torch_dest.device)
    else:
        return torch.tensor(sirf_src.as_array()).to(torch_dest.device)
```

```python
def torch_to_sirf(torch_src, sirf_dest): 
    return sirf_dest.fill(torch_src.detach().cpu().numpy())
```


---
Next we have subclassed `torch.autograd.Function`s these allow for PyTorch to make use of SIRF's acquisition models and objective functions.

Generally they require a combinations of `torch_image`, `torch_measurements`, `sirf_image`, `sirf_measurements`, `sirf_acq_mdl`, and/or `sirf_obj_func`.

In this scope we assume that the dimensionality of arrays is the same across torch/sirf. This therefore uses SIRF's more strict structure, i.e. images of [z,y,x] even if this means there is a singleton dimension.

```python
class _ObjectiveFunctionModule(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            sirf_image_template, 
            sirf_obj_func
            ):

        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        value_np = sirf_obj_func.get_value(sirf_image_template).as_array()
        if torch_image.requires_grad:
            ctx.save_for_backward(torch_image, sirf_image_template, sirf_obj_func)
            return torch.tensor(value_np, requires_grad=True).to(torch_image.device)
        else:
            return torch.tensor(value_np).to(torch_image.device)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        torch_image, sirf_image_template, sirf_obj_func = ctx.saved_tensors
        tmp_grad = sirf_obj.get_gradient(sirf_image_template)
        grad = sirf_to_torch(tmp_grad, torch_image, requires_grad=True)
        return grad_output*grad, None, None, None
```


```python
class _AcquisitionModelForward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_image,
            torch_measurements_template,
            sirf_image_template,
            sirf_acq_mdl
            ):
        
        sirf_image_template = torch_to_sirf(torch_image, sirf_image_template)
        sirf_forward_projected = sirf_acq_mdl.forward(sirf_image_template)
        if torch_image.requires_grad:
            ctx.torch_image = torch_image
            ctx.sirf_forward_projected = sirf_forward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_forward_projected, torch_measurements_template, requires_grad=True)
        else:
            return sirf_to_torch(sirf_forward_projected, torch_measurements_template)

    @staticmethod
    def backward(ctx,
            grad_output
            ):
        sirf_image = ctx.sirf_acq_mdl.backward(torch_to_sirf(grad_output, ctx.sirf_forward_projected))
        grad = sirf_to_torch(sirf_image, ctx.torch_image, requires_grad=True)
        return grad, None, None, None, None
```


```python
class _AcquisitionModelBackward(torch.autograd.Function):
    @staticmethod
    def forward(ctx,
            torch_measurements,
            torch_image_template,
            sirf_measurements_template,
            sirf_acq_mdl
            ):
        
        sirf_measurements_template = torch_to_sirf(torch_measurements, sirf_measurements_template)
        sirf_backward_projected = sirf_acq_mdl.backward(sirf_measurements_template)
        if torch_image_template.requires_grad:
            ctx.torch_measurements = torch_measurements
            ctx.sirf_backward_projected = sirf_backward_projected
            ctx.sirf_acq_mdl = sirf_acq_mdl
            return sirf_to_torch(sirf_backward_projected, torch_image_template, requires_grad=True)
        else:
            return sirf_to_torch(sirf_backward_projected, torch_image_template)

    @staticmethod
    def backward(ctx,
            grad_output
            ):

        sirf_measurements = ctx.sirf_acq_mdl.forward(torch_to_sirf(grad_output, ctx.sirf_backward_projected))
        grad = sirf_to_torch(sirf_measurements, ctx.torch_measurements, requires_grad=True)
        return grad, None, None, None, None
```

---

Next are the classes meant for the user, these subclass `torch.nn.Module`s, and use the aforementioned `Function`s.

There are a few issues. 
- Take for example `torch.nn.conv2d`, there the input/output shape is assumed to be `[N, C, H, W]` or `[C, H, W]`, where the letters correspond to batch (N)umber, (C)hannel, (H)eight, and (W)idth. This needs to integrate into sirf image shape which is just `[z, y, x]`. For `torch.nn.conv3d` this would be `[N, C, D, H, W]` or `[C, D, H, W]`.
- For the measurements this is more complicated, and I am not certain which `torch.nn` operations are even fitting - for example. what does using a trainable convolution on reshaped listmode data really mean...
- The forward operator could be dependent on `N`, i.e. the dataset sample, meaning we'd have to integrate a way of swapping out components of acquisition model (scatter, coil maps, attenuation etc).

Perhaps rather than trying to account for every use case I can make a **simple** `torch.nn.Module` that shows how to manipulate the dimensions to use the `torch.autograd.Function`?

```python
class AcquisitionModelForward(torch.nn.Module):
    def __init__(self, acq_mdl, sirf_image_template, sirf_measurements_template, device = "cpu", ):
        super(AcquisitionModelForward, self).__init__()
        # get the shape of image and measurements
        self.acq_mdl = acq_mdl
        self.sirf_image_shape = sirf_image_template.as_array().shape
        self.sirf_measurements_shape = sirf_measurements_template.as_array().shape
        self.sirf_image_template = sirf_image_template
        self.sirf_measurements_template = sirf_measurements_template
        
        self.torch_measurements_template = torch.tensor(sirf_measurements_template.as_array(), requires_grad=False).to(device)*0

    def forward(self, torch_image):
        # view as torch sometimes doesn't like singleton dimensions
        torch_image = torch_image.view(self.sirf_image_shape)
        return _AcquisitionModelForward.apply(torch_image, self.torch_measurements_template, self.sirf_image_template, self.acq_mdl).squeeze()
```

---
## General notes

What are your thoughts of keeping the `torch.autograd.Function`s follow SIRF's data structure? I think it is more generalisable if we enforce that torch and sirf arrays have the same number of singleton dimensions etc...

Also this means that we can use the `torch.nn.Module`s to do all the data manipulations within the forward pass that is traced. This I think is preferable as messing up something within `torch.autograd.Function` may not cause errors.
