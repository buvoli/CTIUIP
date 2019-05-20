================================================================================

	This directory contains MATLAB Code for:

		Buvoli T. and Tokman M. "Constructing Time Integrators Using 
		Interpolating Polynomials", 2019, Submitted.

	All code is released under the MIT license (See LICENCE).

================================================================================

Directory Structure: 

	ConsantTimestep - Contains abstract classes and implementations for all 
		          integrators.

	Linear Solvers - Contains the abstract classes LinearSolver and two 
			 linear solvers: one for MATLAB mldivide, and one for
			 MATLAB GMRES.

	Nonlinear Solvers - Contains the abstract class NonLinearSolver and
			    a Newton.
	
	Phi - Contains abstract class PhiInitializer and a class for KIOPs.

	Polynomial - Contains two functions for computing coefficients of certain
		     ODE Polynomials.

	Tools - Contains Miscellaneous Helper functions and a class for 
		initializing finite difference matrices.

	Stats - Contains severals Stats Objects.

================================================================================