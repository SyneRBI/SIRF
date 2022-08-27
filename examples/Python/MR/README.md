This folder contains scripts demonstrating the use of SIRF for MR reconstruction and related tasks.

All scripts currently use Gadgetron as the reconstruction engine, however only the demos located in the subfolder named `Gadgetron` employ Gadgetron's gadget chains explicitly. Those outside of this subfolder allow the use of any reconstruction engine that provides the necessary functionality such as `AcquisitionData` object, `preprocess_acquisition_data` function etc.

See the [Python/README](../README.md) for instructions.

<!--
### Fully Sampled Data Reconstructions

The script `fully_sampled_recon.py` shows the reconstruction of fully sampled data without the explicit use of Gadgetron gadgets.

The scripts `Gadgetron/fully_sampled_recon_single_chain.py`, `Gadgetron/fully_sampled_recon_single_chain_short.py` and `Gadgetron/fully_sampled_recon_two_chains.py` show the reconstruction of fully sampled data with explicit use of chains of Gadgetron gadgets. The first of them demonstrates the use of two standard Gadgetron fully sampled reconstruction chains. The second demonstrates the use of Gadgetron's gadgets together with a gadget set provided by SIRF. The last one shows how a Gadgetron reconstruction chain can be split into 2 separate chains - a preprocessing chain and a reconstruction chain - with data filtering applied between these chains.

### GRAPPA Undersampled reconstructions

The script `grappa_basic.py` shows basic reconstruction of data acquired with GRAPPA undersampling without explicit use of Gadgetron gadgets.

The script `Gadgetron/grappa_detail.py` is similar to `grappa_basic.py` but with more annotations in code, output of g-factors and explicit use of gadgets.

The script `grappa_and_steepest_descent.py` shows how to use MR acquisition model to improve the accuracy of GRAPPA reconstruction via the steepest descent iterations.

### Other

The script `acquisition_data.py` shows how to access and manipulate acquisition data.

The script `acquisition_model.py` shows how to generate simulated acquisition data and backproject it into image space (see `grappa_and_steepest_descent.py` demo for an illustration on how this functionality can be employed in the improvement of the accuracy of reconstruction).

The script `coil_sensitivity_maps.py` demonstrates methods for calculating coil sensitivity maps.

This folder contains scripts intended to be run from the command line.
For instance,
```bash
python acquisition_model.py 
```
The demos also have some options. Try
```bash
python acquisition_model.py --help
```
-->