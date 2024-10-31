Code and Data for Ke et al.: Time will tell: the temporal and demographic contexts of plant-soil microbe interactions

**R scripts**
1. ``Figure2_MotivatingAnnualPlant.R`` -- Simulate annual plant competition with temporally varying microbial effects on different demographic processes (Figure 2).
2. ``Figure3_MetaAnalysis.R`` -- Visualized and obtain simple statistics for the meta-analysis (Figure 3).
3. ``FigureBox1_TemporalPatchOccupancy.R`` -- Simulate a patch occupancy model that incorporates the temporal decay of microbial effects (Figure Box 1).
4. ``FigureBox2_AnnualPerennialSensitivity.R`` -- Simulate annual-perennial competition and identify the relative importance of different microbial effects (Figure Box 2).


**Data**
1. ``Data_DemographicReview_CrawfordAndYan.xlsx`` -- Compiled data from Crawford et al. (2019, Ecology Letters) and Yan et al. (2022, PNAS) used in ``Figure3_MetaAnalysis.R``.
2. ``Teste_et_al_2017_Science_Biomass.xlsx`` -- Data from Teste et al. (2017, Science) used to parameterize the model in ``FigureBox1_TemporalPatchOccupancy.R``.
3. ``DecaySimulation.rds`` -- Simulation outputs used to generate Figure Box 1. This is the output from ``FigureBox1_TemporalPatchOccupancy.R`` and is directly provided here as the simulation will take a while.
4. ``ModelFramework_AnnualModel.pdf`` -- Illustration created by BioRender; provided to make Figure 2 fully reproducible.
5. ``ModelFramework_PatchOccupancyModel.pdf`` -- Illustration created by BioRender; provided to make Figure Box 1 fully reproducible.
