Post-flame agglomeration algorithm (PFAL): A hierarchical procedure to produce more realistic soot structures
===

## Introduction
The present code generates fractal structures that resemble soot through aggregation of initial populations of primary particles. We used two stages of discrete element modeling using Langevin dynamics equations based on an earlier works published in [1] and [2]. We calibrated the aggregates between these two stages based on the universal correlation of primary particle size and aggregate size proposed in [3]. The intention for generation of such structures is to explain atypical observations on soot that show in the form of hybridity in primary particle size as reported in [4]. We used the present PFAL package to produce hybrid soot and predicted, for the first time, that those hybrid structures have unique behaviours in terms of projected area and effective density. The results are published in [5]. Examples of real hybrid soot imaged by the author under transmission electron microscope can be found below:

![image](https://github.com/user-attachments/assets/1d6368d9-7990-4174-a418-7c64e2aae5ce)

Discrete element models such as the individual stages of Langeving dynamics modeling used here generally follow a series of simulation steps to produce fractal aggregates:
1. Reading input parameters given by the user related to aggregation phenomena (e.g. number, size and polydispersity of primary particles, volume fraction,...),
2. Random initializing the locations and velocities of primary particles,
3. Calculating the mobility of particles
4. Moving the particles via the interactions of brownian motion, drag and inertia,
5. Imposing periodic conditions on the boundaries (escaping particles return back to the domain from the other side),
6. Checking for particle collisions at every iteration,
7. Attaching the collided particles to form larger aggregates,
8. Calculating the new size and mobility of aggregates to either stop the simulations or move to the next iteration.

To accomplish the above steps, the following packages come with PFAL:
* "+PAR" initializes size, location and velocity of particles and calculates rotation, various measures of size, projected area (via Monte Carlo), densisty-density correlation, etc.
* "+TRANSP" computes transport properties of particles, solves the particle equation of motions and applies the proper boundary conditions,
* "+COL" contains different tools to check collisions and attach the particles to form fractal aggregates,
* "+UTILS" and "+VIS" include a variety of functionalities to visualize and post-process the results,
* "+DEPOT" includes useful functions developed during different phases of this project, but did not end up in final version due to perational reasons.
* "input" contains batch files to initialize the program.

Below are examples of the outputs produced by this program:

![image](https://github.com/user-attachments/assets/6c4902be-6d6a-4ed0-847a-68a19cb2c083)

## References:
1. Suresh, V., & Gopalakrishnan, R. (2021). Tutorial: Langevin Dynamics methods for aerosol particle trajectory simulations and collision rate constant modeling. Journal of Aerosol Science, 155, 105746.
2. Heine, M. C., & Pratsinis, S. E. (2007). Brownian coagulation at high concentration. Langmuir, 23(19), 9882-9890.
3. Olfert, J., & Rogak, S. (2019). Universal relations between soot effective density and primary particle size for common combustion sources. Aerosol Science and Technology, 53(5), 485-492.
4. Baldelli, A., Trivanovic, U., Corbin, J. C., Lobo, P., Gagné, S., Mille, J. W., ... & Rogak, S. (2020). Typical and atypical morphology of non-volatile particles from a diesel and natural gas marine engine. Aerosol and Air Quality Research, 20(4), 730-740.
5. Nikookar, H., Sipkens, T. A., & Rogak, S. N. (2025). Simulating the effect of post-flame agglomeration on the structure of soot. Aerosol Science and Technology, 59(1), 1-15.
