# Fortran

COMPUTATIONAL FLUID DYNAMICS - DIFFUSION AND CONVENCTION TRANSPORTATION IN 1D SPACE

The All the schemes.f95 solves the transportation in one-dimensional space of a property by convection and diffusion means. The domain is first discretized into ten control volumes and then the governing equation is linearized using the Central Differential Scheme, Upwind Differential Scheme and Power Law Differential Scheme. A promp window will ask which scheme to use. The code was written and compiled using the IDE Netbeans. Nevertheless, the code can be compiled using any other platforms.

The Error.f95 calculates the error between the analytical and numerical solution (obtaind using one of the schemes mentioned above). In addition, the code allows the user to explore how numerical error is affected by the number of control volumens. 
