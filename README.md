# BCB330Y1 individual project

-   This repo is for individual scRNA-seq project as part of my BCB330Y1 course
-   Raw data excluded

# Progress

1.  First generate annotated umap with most cluster identified (2025/5/21).\
    Control:\
    ![control](data/processed/Yuzwa/mice_control/mice_control_annotated.png)\
    Treatment:\
    ![treatment](data/processed/Yuzwa/mice_treatment/mice_treatment_annotated.png)

2.  Generate integrated umap, cell identity not yet annotated (2025/6/3).\
    Integrated:\
    ![integrated](data/processed/Yuzwa/mice_integrated/mice_merged_umap.png)\
    And [research proposal](BCB330_Proposal_Jiaqi_Ma.pdf) written last week

3.  Utilize [Clustree](https://github.com/lazappi/clustree) package to find optimal resolution for clustering.\
    And refactored with GO analysis into another helper_function.R to maintain clean workspace.\
    And move all .R file from root to ~/src file (2025/6/6).

4.  Complete annotation by three tools (2025/6/12).\
    My brain:\
    ![Manual](data/processed/Yuzwa/mice_integrated/mice_merged_manual_annotated.png)\
    [Clustifyr](https://github.com/rnabioco/clustifyr):\
    ![Clustifyr](data/processed/Yuzwa/mice_integrated/mice_merged_clustifyr_annotated.png)\
    [SingleR](https://github.com/dviraran/SingleR): (SingleR generate different cluster because it's annotated at single-cell level)\
    ![SingleR](data/processed/Yuzwa/mice_integrated/mice_merged_SingleR_annotated.png)

5. Complete DEGs analysis between control & treatment of Yuzwa's lab WMI dataset (2025/6/25).\
   Manually determined "Top" DEGs:\
   ![Selected_DEGs](data/processed/Yuzwa/mice_integrated/Selected_DEGs.png)\
   Influenced by the recent doubt about p-value == 0.05 initiated from this [Paper](https://doi.org/10.1080/00031305.2016.1154108). I calculate another `p_FC` value, which is a function of p_val and avg_log2FC to help me determine the "Significance" of DEGs rather than thresholding p == 0.05.
   ### `p_FC = f(p_val, avg_log2FC) = (1 - p_val) * (avg_log2FC)`
   
