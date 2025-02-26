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
```

```python
def torch_to_sirf(torch_src, sirf_dest): 
```


---
Next we have subclassed `torch.autograd.Function`s these allow for PyTorch to make use of SIRF's acquisition models and objective functions.

Generally they require a combinations of `torch_image`, `torch_measurements`, `sirf_image`, `sirf_measurements`, `sirf_acq_mdl`, and/or `sirf_obj_func`.

In this scope we assume that the dimensionality of arrays is the same across torch/sirf. This therefore uses SIRF's more strict structure, i.e. images of [z,y,x] even if this means there is a singleton dimension.

```python
class _ObjectiveFunction(torch.autograd.Function):
```


```python
class _OperatorForward(torch.autograd.Function):
```


```python
class _OperatorBackward(torch.autograd.Function):
```

---

Next are the classes meant for the user, these subclass `torch.nn.Module`s, and use the aforementioned `Function`s.

There are a few issues. 
- Take for example `torch.nn.conv2d`, there the input/output shape is assumed to be `[N, C, H, W]` or `[C, H, W]`, where the letters correspond to batch (N)umber, (C)hannel, (H)eight, and (W)idth. This needs to integrate into sirf image shape which is just `[z, y, x]`. For `torch.nn.conv3d` this would be `[N, C, D, H, W]` or `[C, D, H, W]`.
- For the measurements this is more complicated, and I am not certain which `torch.nn` operations are even fitting - for example. what does using a trainable convolution on reshaped listmode data really mean...
- The forward operator could be dependent on `N`, i.e. the dataset sample, meaning we'd have to integrate a way of swapping out components of acquisition model (scatter, coil maps, attenuation etc).

Perhaps rather than trying to account for every use case I can make a **simple** `torch.nn.Module` that shows how to manipulate the dimensions to use the `torch.autograd.Function`?

```python
class ObjectiveFunction(torch.nn.Module):
```

```python
class Operator(torch.nn.Module):
```

---

## Use cases?

- We can do gradient descent with ADAM using either wrapped objective, or acquisition model with torch's pnll function.
- `torch.autograd.gradcheck`, this is very expensive and is a finite difference (central difference) numerical approximation of the Jacobian to compare gradients.
- SIRF tests do due diligence of adjointness etc. We could just check the gradients passed are the same as those in SIRF...?
- Could test <1,Ax> and gradient A^* 1, from both sirf and pytorch.
- Compare objective function value with that computed from wrapped acquistion model and torch's pnll function. 

---
## General notes

What are your thoughts of keeping the `torch.autograd.Function`s follow SIRF's data structure? I think it is more generalisable if we enforce that torch and sirf arrays have the same number of singleton dimensions etc...

Also this means that we can use the `torch.nn.Module`s to do all the data manipulations within the forward pass that is traced. This I think is preferable as messing up something within `torch.autograd.Function` may not cause errors.
