# Pediatric Autoimmunity and Diabetes

In this study we analyze flowCytometry data to characterize the peripheral blood immune profile of children with new-onset T1D. Relatives of T1D patients with 0-1 islets autoantibodies, patients with other autoimmune diseases (celiac disease and autoimmune thyroiditis) and healthy patients were used as controls.
Fresh whole blood samples were stained with five panels of antibodies against 26 surface markers and the intracellular marker FOXP3.

![image](https://user-images.githubusercontent.com/77587985/137467040-095d37dc-5974-4c47-a87c-2a75ed7272a6.png)

For data analysis, unsupervised computational analyses (clustering and visualization) was performed along with statistical analyses to identify immune signatures significantly different between the clinical groups; automated supervised gating analyses were subsequently performed to validate the results.

# FlowCytometry data analysis

you can find a folder containing the main analysis scripts:

* preprocessing

Preprocessing step consists of the import of .fcs files in R environment, the logicle-transformation of fluorescence intensities and compensation of spillover coeffiencts. Spillover coefficients were estimated for each fluorophore using single-color controls with autospill algorithm [Carlos P. Roca, Oliver T. Burton, Teresa Prezzemolo, Carly E. Whyte, Richard Halpert, Łukasz Kreft, James Collier, Alexander Botzki, Josef Spidlen, Stéphanie Humblet-Baron and Adrian Liston. AutoSpill: a method for calculating spillover coefficients in high-parameter flow cytometry].
To adjust for between-sample and batch variation, all markers satisfying the rules for normalization, were subject to the landmark alignment procedure. Finally, cells have been gated to extract the cells populations used as input for flowSOM algorithm. Gatings of cells were performed using OpenCyto R package.

* unsupervised_analysis

The unsupervised analysis has been done with flowSOM algoritm included in the CyTOF / CATALYST pipeline 
[Crowell, H. L. et al. An R-based reproducible and user-friendly preprocessing pipeline for CyTOF data. F1000Res. 9, 1263 (2020)]

* supervised_gating_strategy

Supervised analysis was performed through automated gating R package OpenCyto 
[Finak G, Frelinger J, Jiang W, et al. OpenCyto: an open source infrastructure for scalable, robust, reproducible, and automated, end-to-end flow cytometry data analysis. PLoS Comput Biol. 2014;10(8):e1003806. Published 2014 Aug 28. doi:10.1371/journal.pcbi.1003806]


# Statistics 

In the downstream_analysis folder you can find the scripts for producing plots and statistics.


