Welcome!

This is **GIZMO**: a flexible, multi-method multi-physics code. The code solves the fluid using Lagrangian mesh-free finite-volume Godunov methods (or SPH, or fixed-grid Eulerian methods), and self-gravity with fast hybrid PM-Tree methods and fully-adaptive resolution. Other physics include: magnetic fields (ideal and non-ideal), radiation-hydrodynamics, anisotropic conduction and viscosity, sub-grid turbulent diffusion, radiative cooling, cosmological integration, sink particles, dust-gas mixtures, cosmic rays, degenerate equations of state, galaxy/star/black hole formation and feedback, self-interacting and scalar-field dark matter, on-the-fly structure finding, and more. 

See the [User Guide](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html) for an up-to-date physics list, or [the original code website](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO.html) for examples demonstrating the advantages of the new methods, different types of **GIZMO** simulations, and its massively-parallel scalings.

The code is descended from P-SPH, itself descended from GADGET-3 (so a huge debt owes to the work of Volker Springel), and many of the GADGET conventions remain (for compatibility with GADGET outputs and codes). See the source code for appropriate attribution of the code elements. 

Both the [public code](https://bitbucket.org/phopkins/gizmo-public) and the [development (private) code](https://bitbucket.org/phopkins/gizmo) are hosted on Bitbucket version-controlled repositories, but for code issues, feature requests, or bugs, please post to the [GIZMO Google Group](https://groups.google.com/d/forum/gizmo-code).

Basic Rules: 

1. The reference for code methods, setup, use and citation policy is the User Guide, available as `gizmo_documentation.html` in the `scripts` folder as part of this repository, or via download [here](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html). Read it! The code is extensively documented. Important things you need to know before running are there. Most questions are already answered.  

2. Access to the development (private) code does **not** imply permission to use any modules identified as proprietary either in the User Guide or `Template_Config.sh` file (or elsewhere). The development code can only be used or distributed with explicit permission from the code authors. Many of the non-public modules are proprietary and developed by students for on-going research; it is not acceptable to use or share these routines without first obtaining the explicit permission of both the lead code author and author(s) of the relevant routines. If in doubt, ask. Anyone violating these terms will have code access immediately revoked.

The public version of the code is free software, distributed under the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html). You may freely distribute and copy the public code. You may also modify it as you wish, and distribute these modified versions as long as you indicate prominently any changes you made in the original code, and as long as you leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. Note that the authors retain their copyright on the code. 

If you use any version of the code, please reference the code paper, [Hopkins 2015, MNRAS, 450, 53](http://arxiv.org/abs/1409.7395); you should also reference Volker Springel's GADGET paper (Springel, 2005, MNRAS, 364, 1105) for the domain decomposition and N-body algorithms. Appropriate citations for specific modules are described in the Users Guide.
