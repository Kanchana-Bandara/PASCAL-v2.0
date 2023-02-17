## Pan-Arctic Behavioural and Life-history Simulator, PASCAL (v2.00)
#### Species-specific (ensembel) behavioural and life-history simulation model for north Atlantic and Arctic Calanus spp. (1st updgrade to v.1.00)
This repository contains the second version of PASCAL, initiated in 2019 & completed in 2020. PASCAL v.2.00 contains a species-specific, high-resolution behavioural and life-history simulation model for _C. finmarchicus_, _C. glacialis_ and _C. hyperboreus_. 

The outputs from this model have already been analyzed and published once at:

_Bandara, K., Varpe, Ø., Ji, R., & Eiane, K. (2019). Artificial evolution of behavioral and life history strategies of high-latitude copepods in response to bottom-up and top-down selection pressures. Progress in Oceanography, 173, 134-164._

_Bandara, k. (under review) The simulated overlap of body sizes between C. finmarchicus and C. glacialis under bottom-up and top-down selection pressures and its taxonomy implication. PLOS Computational Biology (PCB202212987)._

##### Contents
This repository contains two main files: (i) the PASCAL v2.00 hull files for the 3 species and (ii) PASCAL v2.00 simulation drive for all species, (iii) model environment generation files and (iv). two functions to generate body size trajectories for simulated evolution. The hull file is the main file from which all other functions (model environment, size trajectories, stage-specific behavioural & life-history simulations) are called. Similar to v.1.00, the additional functions file contains three (3) small functions to be loaded besides the two main files.

##### Warnings
The environment files (2D arrays) are needed to execute the model. These files are not included here. Check the corresponding Zenodo repository for example files.
The HPC framework has been redarcted due to an ongoing publication. We will add this as soon as we can in a separate repository (this README file will be updated when the HPC framework is made Open Access)

###### Runtime
The model is set-up to run for 400 generations (calendar years). With the HPC framework built-in, it takes about 42-48 hours to complete execution in a CORE i9 7920X 24-Thread, 2.9-3.5 GHz gaming rig. If you have the PushoverR app, the model will send an update to your phone once the simulation has been completed. For single threaded operations (no HPC), please dial down the no. of individuals in the model. However, the impact of the simulation population size on the emergent ecological properties has not been assessed. 

##### Lisence
CC BY 4.0
