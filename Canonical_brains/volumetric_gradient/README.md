# Volumetric gradient maps

We (Jae-Joong Lee and Suhwan Gim) estimated volumetric gradient maps with similar parameters to the original paper [(Margulies et al., 2016)](https://www.pnas.org/doi/10.1073/pnas.1608282113). In the beginning, we used the second version (N=56). However, we recommend you use the first version (N=986) for more generalizability. 


### A volumetric functional connectivity gradient map, based on HCP S1200 dataset (*N* = 986)
   - Default mask used in RF-ANTs (>50% overlap across all Genomics Superstruct Project [GSP] participants; see original RF-ANTS paper: [Wu et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29770530/))
   - Custom script (DM to J.J. if you need), 90% sparsity, diffusion embedding algorithm
   - Brain maps (from gradient 1 to gradient 10): https://github.com/cocoanlab/cocoanCORE/blob/master/Canonical_brains/Volumetric_hcp_gradients_GSP_90_DE_wholebrain_2mm.nii


### A volumetric functional connectivity gradient map using our own dataset (*N* = 56)

   - Average 56 functional connectivities of a resting-state  fMRI (3mm iso voxel, 810 volumes, around 6 min and 10 secs, TR = 0.46 secs, parcellation using brainnetome atlas)
   - BrainSpaceToolbox in MATLAB (Normalized angle, 90% sparsity, diffusion embedding)
   - Brain maps (from gradient 1 to gradient 10): https://github.com/cocoanlab/cocoanCORE/blob/master/Canonical_brains/ten_prinicpal_gradients_volumn_DE.nii
   - Data from [Gim et al., 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002910).



### Contacts 


Jae-Joong Lee (jaejoonglee92@gmail.com) & Suhwan Gim (suhwan.gim.psych@gmail.com)

___
### Example (see Figures 4 and 5 in Gim et al., 2024)
You can calculate overlapped proportions between your map (thresholded) with a principal gradient map using this function [principal_gradients_plot.m](https://github.com/cocoanlab/cocoanCORE/blob/master/Visualization/principal_gradients_plot.m) in [cocoanCORE](https://github.com/cocoanlab/cocoanCORE/) 
(You may need other toolboxes such as [canlabCORE](https://github.com/canlab/CanlabCore)) 

### References

Gim, S., Hong, S. J., Reynolds Losin, E. A., & Woo, C. W. (2024). Spatiotemporal integration of contextual and sensory information within the cortical hierarchy in human pain experience. PLoS biology, 22(11), e3002910

Margulies, D. S., Ghosh, S. S., Goulas, A., Falkiewicz, M., Huntenburg, J. M., Langs, G., ... & Smallwood, J. (2016). Situating the default-mode network along a principal gradient of macroscale cortical organization. Proceedings of the National Academy of Sciences, 113(44), 12574-12579.

Wu, J., Ngo, G. H., Greve, D., Li, J., He, T., Fischl, B., ... & Yeo, B. T. (2018). Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems. Human brain mapping, 39(9), 3793-3808.
