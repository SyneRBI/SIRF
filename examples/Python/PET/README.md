This folder contains scripts demonstrating the use of SIRF for PET reconstruction and related tasks.

All scripts currently use STIR as the reconstruction engine but actually allow the use of any reconstruction engine that provides the necessary functionality such as `AcquisitionData` object, `make_Poisson_loglikelihood` function etc.

<!--
This folder contains scripts intended to be run from the command line.
For instance,
```bash
python acquisition_model.py 
```
The demos also have some options. Try
```bash
python acquisition_model.py --help
```

Check the `interactive` folder for demos that you run from an IDE.
-->