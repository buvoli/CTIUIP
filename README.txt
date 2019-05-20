================================================================================

	This directory contains MATLAB Code for:

		Buvoli T. and Tokman M. "Constructing Time Integrators Using 
		Interpolating Polynomials", 2019, Submitted.

	All code is released under the MIT license (See LICENCE).

================================================================================

Directory Structure: 

	Integrators - Contains abstract classes and implementations for all 
		      integrators.

	Problems - Contains abstract classes for problem, and classes for the
		   2D ADR equation and the 1D Burgers' equation.

	Tools - Contains Miscellaneous Helper functions and a class for 
		initializing finite difference matrices.

================================================================================

Example Usage:
	
	Run main_example to compare BDF, AM, BBDF and BAM on either Burgers' or
	2D ADR. 