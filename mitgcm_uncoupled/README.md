# README

This code and input folder are the configuration for MITgcm to be used in a first test for coupling.

The configuration is based on 

Jordan, J. R., Holland, P. R., Goldberg, D., Snow, K., Arthern, R., Campin, J.-M., … Jenkins, A. (2018). Ocean-forced ice-shelf thinning in a synchronously coupled ice-ocean model. Journal of Geophysical Research: Oceans, 123, 864– 882. https://doi.org/10.1002/2017JC013251

Nakayama, Y., Hirata, T., Goldberg, D., & Greene, C. A. (2022). What determines the shape of a Pine-Island-like ice shelf? Geophysical Research Letters, 49, e2022GL101272. https://doi.org/10.1029/2022GL101272

And consists of an ice shelf with an upstream boundary condition flowing into a rectangular cavity. In this configuration streamice is not used. 

To produce input binary files for the experiment, run the following files in this order:
`rdmds_init.m`
`gendata.m`
