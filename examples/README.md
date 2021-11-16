# Examples

In order to use these codes, please use the Docker container file ``../Dockerfile``) which can be build according to the README (https://github.com/LFPy/LFPykernels#docker).
This cross-platform solution provides an environment including all dependencies for the main example notebook ``LIF_net_forward_model_predictions.ipynb``

## Usage

Assuming the Docker container has been built according to this project's README, the example notebook(s) may be executed by issuing in the terminal:

    cd <path-to-LFPykernels>
    docker run --mount type=bind,source="$(pwd)",target=/opt/data -it -p 5000:5000 lfpykernels
    root@b6ff6a542d8a:/# cd /opt/data/examples
    root@b6ff6a542d8a:/# jupyter-notebook --ip 0.0.0.0 --port=5000 --no-browser --allow-root

Then connect to the server with URL similar to http://127.0.0.1:5000/?token=6c26f9a5a9c18f52c31a572ba3bda255f278a40a91297a55 using a browser on the host computer. There provided example notebook(s) may be run and edited in an interactive manner.

## File list

- `BallAndSticks_E.hoc`
  Excitatory cell morphology file
- `BallAndSticks_I.hoc`
  Inhibitory cell morphology file
- `BallAndSticksTemplate.hoc`
  Cell template specification file needed by `LFPy.NetworkCell` loading the different morphologies
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
- `README_example.ipynb`
  Reference implementation of example code in main README file. 
