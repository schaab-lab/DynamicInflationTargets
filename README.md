# DynamicInflationTargets

This repository accompanies our paper:

Clayton, C. and A. Schaab. A Theory of Dynamic Inflation Targets.Conditionally accepted at American Economic Review, 2024.

The script main.m computes the impulse response functions for our applications in Sections 4, 5 and Appendix B.2. In particular, it generates Figures 1 (Section 4.1), 2 (Section 4.2) and 3 (Section 5.3), as well as Figure 4 (Appendix B.2). These Figures are printed to the folder /output.

Our script uses an axis break routine for Figure 3:
Chris McComb (2024). Aesethetic Axis Breaks (https://github.com/cmccomb/break_axis), GitHub.

# Script and functions:
- main.m is the master script that generates all Figures
- define_parameters.m is a function that defines all relevant economic and grid parameters
- figure_format.m is a script that initializes Figure format options
- figure_irf_panel.m is a function that formats Figures
- break_axis.m is a function developed by Chris McComb for an aesthetic axis break (Figure 3)
