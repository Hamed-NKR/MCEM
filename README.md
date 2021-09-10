Monodisperse clustering external mixing (MCEM) program: An algorithm for more realistic aggregation modeling of soot
===

The present code aims to generate fractal soot aggregates from clustering of primary particle populations. The aggregation phenomenon studied here is in the diffusion-limited cluster-cluster aggregation (DLCA) regime. The numerical algorithm used consists of solving transport equation for the particles in a Langevin Dynamics (LD) framewrok [1] and simulating their growth upon collision with each other. The primary particles considered are nascent soot in typical flame conditions having a polydisperse size distribution. The long-term goal for this project is to model the aggregation in a new way which assumes a hierarchical clustering suggested by External Mixing Hypothesis [2]. This will results in aggeragtes which are externally polydisperse (in terms of primary particle size) while having an internally monodisperse morphology.

LD-DLCA algorithms such as MCEM typically follow the below steps to generate agglomerated soot [3]:
1. Reading the aggregation input parameters,
2. Random initialization of the monomers' locations and velocities,
3. Moving the particles based on brownian motions and drag,
4. Imposing periodic conditions on the boundaries,
5. Checking for particle collisions over time,
6. Attaching the collided particles as larger clusters,
7. Maintianing the volume fraction by enlarging the main computational domain upon clusteraions,
8. Exporting and plotting the aggregation data.

To accomplish the above steps the following packages come with MCEM:
* "+PAR" initializes different properties of particles,
* "+TRANSP" computes important transport properties of fluid and particles, solves the particle equation of motions and applies the proper boundary conditions,
* "+COL" contains different tools to check collisions and cluster the particles,
* "+UTILS" implements different tools to visualize and post-process the results.
* "+DEPOT" includes useful tools generated through development of the program, but not employed in the final version for different reasons.

In addition to these, the "input" and "output" folders, respectively, contain the initialization parameters and final results of the program. To run the program, two scripts are included in the main directory named as "main_test" and "main_fast". The former one is written to monitor and troubleshoot the outputs while the latter aims to run the code in the lightest form possible to get the results at the lowest computational cost. The program can also be run in two different modes of data storage:
* via a simpler global structure having a concatinated form of aggregates information. This is expected to be faster as a result of to lower amounts of looping, but trickier to develop due to complexity,
* with a class of aggregates containing different global and internal information of each aggregate object within them. This is more intuitive and easier to deal with, but requires more computational demand.

# References:
1. Suresh, V., & Gopalakrishnan, R. (2021). Tutorial: Langevin Dynamics methods for aerosol particle trajectory simulations and collision rate constant modeling. Journal of Aerosol Science, 155, 105746.
2. Olfert, J., & Rogak, S. (2019). Universal relations between soot effective density and primary particle size for common combustion sources. Aerosol Science and Technology, 53(5), 485-492.
3. Heine, M. C., & Pratsinis, S. E. (2007). Brownian coagulation at high concentration. Langmuir, 23(19), 9882-9890.
