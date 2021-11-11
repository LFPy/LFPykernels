# Examples

In order to use these codes, please use the Docker container (`../Dockerfile`) which can be build according to the README (https://github.com/LFPy/LFPykernels#docker).


## File list

- `BallAndSticks_E.hoc`
  Excitatory cell morphology
- `BallAndSticks_I.hoc`
  Inhibitory cell morphology
- `BallAndSticksTemplate.hoc`
  Cell template specification needed by `LFPy.NetworkCell` loading the different morphologies
- `FIR_filter.nestml`
  NESTML (https://github.com/nest/nestml) file specifying a Finite Impulse Response (FIR) filter node in NEST
- `example_network_methods.py`
  Various methods required by network simulations
- `example_network_parameters.py`
  Parameters for recurrent multicompartment network model set up using LFPy. Here used to derive parameters for kernel predictions
- `plotting.py`
  Some shared plotting methods.
- `mod/*.mod`
  NMODL language files describing ion-channel and synapse dynamics
- `LIF_net_forward_model_predictions.ipynb`
  Jupyter notebook implementing a spiking point-neuron network model with forward-model based predictions using the `LFPykernels` python package
