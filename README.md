# Kv-kinetic-models
Kinetic models of voltage-gated potassium channels.
For every channel different models were fit to experimental data kindly provided
by members of the Blue Brain Project. The models were fit directly to the experimental data
using the data2dynamics software (https://github.com/Data2Dynamics/d2d).

The models are organized with one directory for each channel and a subdirectory
for each model. The subdirectories contain the model, a description of the model,
a list of the model parameters, a summary of the behavior of the model under different
experimental conditions and a NEURON hoc file to test the model. 

For each channel there is also a subdirectory containing an example trace for propagation
of an action potential in a standard model of a neuronal cell, but using the given channel model
for the potassic current, i.e. we replace the normal potassic channel of the model (a Hodgkin-Huxley model)
with one of our models for a specific potassium channel.
An example of how to interpret these figures is given here.
We look at the file "hbp-00009_Kv2.1_AP_propagation_13States_temperature2_1.00-10.0_axon.png"
In the filename we have the potassic channel used (KV2.1), the model used (13States_temperature2)
and then 2 numbers. The first one is always 1.00 here, meaning that the considered potassium 
channel model is the only potassic current present in the simulated cell. The second represents the
density of channels used, here 10 times the density of the HH potassium channel model it replaces. The
higher density used is to compensate for some models that yield lower potassium conductances (for example only
a fraction of the channels open).
The figure itself has two panels. In the top panel, we see the voltage traces:
* We start at equilibrium membrane potential
and after 100ms, a continuous current is injected in the soma during a few ms. The transmembrane voltage rises in the soma (red curve) and eventually an AP is fired and propagates down the axon (blue curve is at the beginning of the axon and green curve at the end of the axon).
* In the second panel we show the open state probabilities for the different channels, in the soma (red), at the entrance of the axon (blue) and at the end of the axon (green). The dashed lines show the open probabilities of the sodium channels, a Hodgkin-Huxley model, while the continuous lines show the open probabilities for the potassium channel. For comparison we show in red the open probability of the HH model in the soma, while blue and green show the open probability of the kinetic model of Kv (here the 13state model of Kv2.1). Because this kinetic model has somewhat lower open probabilities, we show 7 times the open probability here (see legend 7.0*kv axon0 and 7.0*kv axon1). In this example it is obvious that the channel closes very slowly compared to the HH model. This leads to a longer Hyperpolarized phase after the AP.