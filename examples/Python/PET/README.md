This folder contains scripts demonstrating the use of SIRF for PET reconstruction and related tasks.

All scripts currently use STIR as the reconstruction engine but actually allow the use of any reconstruction engine that provides the necessary functionality such as `AcquisitionData` object, `make_Poisson_loglikelihood` function etc.

See the [Python/README](../README.md) for instructions.


| demo script  | purpose |
|--|--|
|`acquisition_data.py`  | illustrates acquisition data handling|
|`acquisition_data_from_scanner_info.py`  | shows how to create acquisition data from scanner parameters and axial compression etc.|
|`acquisition_model.py` | illustrates the use of PET acquisition model: creates an image, projects it to simulate acquisition data and backprojects simulated data|
|`acquisition_sensitivity_from_attenuation.py`  | illustrates acquisition sensitivity model representing attenuation factor|
|`acquisition_sensitivity_from_bin_efficiencies.py`  | illustrates acquisition sensitivity model representing bin efficiencies factor|
|`acquisition_sensitivity_from_ecat8.py`  | illustrates acquisition sensitivity model based on ECAT8 bin normalisation|
|`fbp2d_reconstruction.py`  | illustrates reconstruction by two-dimensional filtered backprojection|
|`get_multiplicative_sinogram.py`  | illustrates obtaining  multiplicative sinograms from normalisation and/or attenuation|
|`hkem_reconstruction.py`  | illustrates Hybrid Kernelized Expectation Maximization reconstruction|
| input_output.py | illustrates reading/writing acquisition data and images from/to files|
|`listmode_to_sinograms.py`  | illustrates the conversion of listmode data to sinograms|
|`osem_reconstruction.py`  | illustrates reconstruction by Ordered Subsets (OS) version of the One Step Late (OSL) algorithm of Green et al for Maximum a Posteriori (MAP) maximisation|
|`osem_reconstruction_gpu.py` | implements the above OSMAPOSL reconstruction on GPU|
|`osl_reconstruction.py` | illustrates reconstruction by OSMAPOSL algorithm with penalty term|
|`osl_reconstruction.py` | illustrates reconstruction by Ordered Subsets Separable Paraboloidal Surrogate (OSSPS) algorithm|
|`randoms_from_listmode.py`  | illustrates estimation of randoms from listmode data|
|`reconstruct_from_listmode.py`  | illustrates the conversion of listmode data to sinograms with subsequent reconstruction|
|`scatter_estimation.py`  | illustrates estimation of scatter|
|`scatter_simulation.py`  | illustrates simulation of scatter|
|`steepest_ascent.py` | applies few steps of steepest ascent for the maximization of Poisson log-likelihood objective function using subset gradients|
|`STIR_acquisition_model_using_raytracing` | illustrates setting advanced parameters in STIR acquisition model that uses raytracing|
|`osl_reconstruction.py` | Python implementation of OSMAPOSL algorithm with penalty term|
