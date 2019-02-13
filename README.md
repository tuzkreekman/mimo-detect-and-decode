# MIMO Decoding

Some very basic code for working with MIMO

## SimZFMMLReal.m

Simulates a few MIMO decoders over real channels with MPAM constellations.
Includes full ML decoding but this won't be called if n is set higher than
4 as it is expensive.

## SimZFMLComplex.m

Same as above but for complex channels and MQAM.  Not much is different here,
you mainly just split of real vs. complex at the receiver when performing
detection

## sphereDecode.m

A basic sphere decoder for square channels

## estimate_channel.m

Performs least squares channel estimation --- this is what the paper 
``Deep Learning for Joint MIMO Detection and Channel Decoding'' claims for their 
baseline comparisions.

## jakes_model.m

Implements a basic SISO time-varying Jakes' scattering model. 

## cartprod.m

This just returns the cartesian product of two vectors.  I use this to generate
QAM constellations
