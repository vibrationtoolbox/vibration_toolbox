% Student Edition of the Vibration Toolbox.
% Version 5.0 1-13-00
% For use with MATLAB 4.x-5.x
% Updates are reflected in vtbud.m
% Also see vtb
%
% General
%   Contents   - contents of vtoolbox
%   vtbud      - "help vtbud" shows updates to vtoolbox since 1/1/98
%
% VTB1: Free response of single degree of freedom (SDOF) systems.
%   vtb1_1     - plots the free response of a single-degree of freedom
%                system.
%   vtb1_2     - numerically integrates to find the unforced response 
%                of a system using Euler's method.  The system can
%                be in state space or LSOM form.
%   vtb1_3     - performs fourth Runge-Kutta integration on a set of
%                equations described in a specified file.  The 
%                equations need not be linear.
%   vtb1_4     - numerically integrates to find the forced response 
%                of a system using Euler's method.  The system can
%                be in state space or LSOM form.
%   vtb1_5     - plots the free decay of a single degree of freedom
%                system with viscous, coulomb, or air damping.
%
% VTB2: Response to harmonic excitation, SDOF.
%   vtb2_1     - Harmonic response of an undamped SDOF system.
%   vtb2_2     - Harmonic response of an underdamped SDOF system.
%   vtb2_3     - FRF (transfer function) of a SDOF system.
%   vtb2_4     - Displacement and force transmissibility ratios.
%   vtb2_5     - Rotating imbalance.
%   vtb2_6     - Free decay with viscous, coulomb, or air damping.
%
% VTB3: General forced response.
%   vtb3_1     - Impulse response of a SDOF system.
%   vtb3_2     - Step response of a SDOF system.
%   vtb3_3     - Fourier Series approximation to a triangle wave.
%   vtb3_4     - Response spectrum for a SDOF system.
%
% VTB4: Multiple-degree-of-freedom systems.
%   vtb4_1     - Frequencies and eigenvectors/undamped system.
%   vtb4_2     - Free response of an undamped system.
%   vtb4_3     - Solves for the natural frequencies, damping ratios, 
%                and mode shapes of a linear second order matrix (LSOM)
%                form system.  (M,C,K)
%
% VTB5: Design for vibration suppression.
%   vtb5_1     - Transmissibility ratio of a SDOF system.
%   vtb5_2     - Base excitation force transmissibility.
%   vtb5_3     - Motion of primary mass for absorber design.
%   vtb5_4     - Mass ratio vs. normalized frequency for absorber.
%   vtb5_5     - Damped vibration absorber design.
%   vtb5_6     - Surface plot of vtb5_5.
%
% VTB6: Distributed-parameter systems.
%   vtb6_1     - Natural frequencies and mode shapes of a uniform bar.
%   vtb6_2     - Torsional frequencies and mode shapes of a uniform bar.
%   vtb6_3     - Natural frequencies of an Euler-Bernoulli beam.
%   vtb6_4     - Comparison of natural frequencies for beams.
%
% VTB7: Vibration testing and experimental modal analysis.
%   vtb7_1     - Power spectral density of data.
%   vtb7_2     - Discrete frequency response function (H(iw)).
%    [vtb7_2ex   - example data file for vtb7_2]
%   vtb7_3_*   - Experimental data files.
%   vtb7_4     - SDOF curve fitting.
%   vtb7_5     - Transfer function calculation from MDOF matrices.
%    
% VTB8: The finite element method.
%   vtb8read   - FEM Manual in text format.
%   vtb8_1     - FEM pre-processor.
%   vtb8_2     - FEM.
%   vtb8_3     - FEM post-processor.
%    [e8_2_1.*   - example 8.1.2, .exa is a diary file ]
%    [p8_3_10.*  - problem 8.3 solved with 10 elements ]
%    [p8_12.m    - problem 8.12                        ]
%    [truss1.*   - truss example                       ]
%    [truss2.*   - truss examples                      ]
%    [vtb8_e1.*  - example problem                     ]
%
% VTB9: Computational considerations.
%   vtb9_1     - Natural frequencies, damping ratios, mode shapes.
%   vtb9_2     - Unforced response of a system. Euler's method.
%   vtb9_3     - Fourth order Runge-Kutta integration.
%    [vtb9_3ex   - example for use with vtb9_3 ]
%   vtb9_4     - Forced response of a system: Euler's method.
%
% VTB10: Nonlinear vibration.
%   vtb10_1    - Plots phase planes.
%    [vtb10ex    - Example function file for vtb10_1]
%
% Type help vtb# for help on toolbox programs for chapter #.
% Type help "filename" for specific help on that file.
% Type help vtbud for a list of updates since 1/1/98.
%
%  Please connect to the Engineering Vibration Toolbox home
%  page at
%  http://www.cs.wright.edu/people/faculty/jslater/vtoolbox/vtoolbox.html
%  for the latest information.


