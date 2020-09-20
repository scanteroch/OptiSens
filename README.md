# OptiSens - Convex optimization of sensor and actuator placement for ultrasonic guided-wave based SHM

`OptiSens` is comprised of two computational platforms (i.e. one in Python and the other one in Matlab) to provide optimal sensor and actuator configurations for structural health monitoring (SHM) applications using ultrasonic guided-waves. The Python version is oriented to developers and researchers who may want to continue developing the source code, the Matlab version is devoted to end-users, such as MSc/PhD students and practitioners.

The `OptiSens` software formulates a convex entropy-based objective function which is minimized, thus minimizing the uncertainty while maximizing the expected accuracy of the resulting monitoring system in localizing damage. The platforms are specialized for two types of different materials, namely isotropic and composite (anisotropic) materials.

## Example

The `OptiSens` software contains several examples when running the default values provided within the scripts for the two types of materials that is oriented: isotropic materials and anisotropic (composite) materials. The output of the isotropic material, considering aluminum properties, is illustrated as follows:

<img src="https://github.com/scanteroch/OptiSens/blob/master/Python/Metals/res/OptSol.png" alt="drawing" width="400"/><img src="https://github.com/scanteroch/OptiSens/blob/master/Matlab/Metals/res/Opt_Sol.png" alt="drawing" width="430"/>

Where the left panel represents the output given by the Python-based platform, while the right panel shows the results provided by the Matlab-based platform. Note that the blue (left-pointing) triangles represents the optimal location for sensors, while the orange (right-pointing) triangles depicts the location of the actuators. Besides, the larger the size of the markers, the best positions for either sensors or actuators. The quantitative and more accurate information is written in the text files, as follows:

<img src="https://github.com/scanteroch/OptiSens/blob/master/Python/Metals/res/Opt_Sol_txt.png" alt="drawing" width="300"/>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://github.com/scanteroch/OptiSens/blob/master/Matlab/Metals/res/Opt_Sol_txt.png" alt="drawing" width="300"/>

In the tables shown above, the most relevant locations for sensors and actuators are listed along with their coordinates and value of decision variable.

## References

The `OptiSens` software is based on the following references:

> Sergio Cantero-Chinchilla, James L. Beck, Manuel Chiachío, Juan Chiachío, Dimitrios Chronopoulos, and Arthur Jones. Optimal sensor and actuator placement for structural health monitoring via an efficient convex cost-benefit optimization. [Mechanical Systems and Signal Processing 144 (2020) 106901](https://doi.org/10.1016/j.ymssp.2020.106901).

> Bhattacharyya, P, Beck, J. Exploiting convexification for Bayesian optimal sensor placement by maximization of mutual information. [Struct Control Health Monit. 2020; 27:e2605](https://doi.org/10.1002/stc.2605).



## Acknowledgements


This work is part of the SAFE-FLY project that has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 721455.
