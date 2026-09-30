# Figure 1C–D Public Demonstration Code

Entry point: `main_Figure1CD_PublicDemo.m`. Set this folder as the MATLAB current working directory, then run:

```matlab
run('main_Figure1CD_PublicDemo.m')
```

The parameters at the beginning of the script specify the wavelength, SLM dimensions, pixel pitch, focal length, magnification, input beam waist, target diameter, sample-plane scale, super-Gaussian order, and iteration count. The main script includes the input beam, 15 μm target, initial phase, SPOT iterations, conventional GS calculation, performance metrics, and plotting routines.

With `useArchivedGSReference = true`, the script loads the GS result used in the manuscript. Set this option to `false` to recalculate GS using a random initial phase; the resulting metrics may vary with the random state. No fixed random seed is specified, and loading the reference result does not depend on a random seed.
