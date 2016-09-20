# Comparison_plot

This respo contains scripts to plot different statistics and compare model performances using the database of isotropic turbulence.

## Dependencies

TODO: Describe the list of Dependencies

## Usage

TODO: Write usage instructions

## Contents

1. `downsample4_Fourier.m`: Load 37 fields (of original cube 364^3), then use ideal Fourier filter to downsample to 96^3 data fields, then save to final data of size 37 x 96 x 96 x 96
2. `energyloss_variousratios.m`: Estimate the energy losses of different subsampling ratios in either 2D space (spanwise-vertical) and time (streamwise)

## History

25/07: Added `downsample4_Fourier.m` and `energyloss_variousratios.m`
