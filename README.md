# SOFC-H2-parameter-study
COMSOL models for the study of the effects of reference exchange current density and electrode ionic conductivity on solid oxide fuel cell (SOFC) operation, as shown in our paper (DOI: 10.1016/j.jpowsour.2024.234493).

<p align="center">
<img src="https://github.com/green-ammonia-institute/SOFC-H2-parameter-study/blob/main/Images/BV-MS%20-%20xH2.png" width=50% height=50%>
</p>

## Multiphysics model description

The models describe a three-dimensional stationary and isothermal electrolyte-supported SOFC using pure hydrogen as fuel. Both Fick's law (FL) and Maxwell-Stefan (MS) mass transport models are included. A detailed explanation of the physics involved can be found in the article.

## Folder structure

Folder [Parameter estimation/](Parameter%20estimation/) contains a model that searches for preliminary values of anodic and cathodic reference exchange current densities. [Mesh analysis/](Mesh%20analysis/) contains the results used for mesh analysis. The files in [i0ref analysis/](i0ref%20analysis/) focus on sensitivity analysis of anodic reference exchange current density for FL and MS models. Likewise, [sigma analysis/](sigma%20analysus/) contains models that study sensitivity with respect to electrode ionic conductivity. In every folder there is a Results subfolder that provides convergence, current and power data for each simulation in `.csv` format.

## Running the software

The models have been implemented in COMSOL Multiphysics® 6.2 with the Fuel Cell & Electrolyzer module, using the following interfaces: Thermodynamics, Chemistry, Secondary Current Distribution, Transport of Concentrated Species, Free and Porous Media Flow, and Reacting Flow. It may be possible to run the models with other COMSOL installations if they provide equivalent functionality, but this has not been tested.

MATLAB files should be placed in the same folder as the `.csv` files in their respective Results subfolder. They calculate relevant quantities and generate vectorized plots.

## Acknowledgements

This work was funded by the Millennium Institute on Green Ammonia as Energy Vector - MIGA (ICN2021\_023) supported by the Millennium Scientific Initiative by the Agencia Nacional de Investigación y Desarrollo (ANID).
