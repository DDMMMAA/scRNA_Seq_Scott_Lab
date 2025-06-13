# BCB330Y1 individual project

-   This repo is for individual scRNA-seq project as part of my BCB330Y1 course
-   Raw data excluded

# Progress

1.  First generate annotated umap with most cluster identified (2025/5/21).\
    Control:\
    ![control](data/processed/mice_control/mice_control_annotated.png)\
    Treatment:\
    ![treatment](data/processed/mice_treatment/mice_treatment_annotated.png)

2.  Generate integrated umap, cell identity not yet annotated (2025/6/3).\
    Integrated:\
    ![integrated](data/processed/mice_integrated/mice_merged_umap.png)\
    And [research proposal](BCB330_Proposal_Jiaqi_Ma.pdf) written last week

3.  Utilize [Clustree](https://github.com/lazappi/clustree) package to find optimal resolution for clustering.\
    And refactored with GO analysis into another helper_function.R to maintain clean workspace.\
    And move all .R file from root to ~/src file (2025/6/6).

4.  Complete annotation by three tools (2025/6/12).\
    My brain:\
    ![Manual](data/processed/mice_integrated/mice_merged_manual_annotated.png)\
    [Clustifyr](https://github.com/rnabioco/clustifyr):\
    ![Clustifyr](data/processed/mice_integrated/mice_merged_clustifyr_annotated.png)\
    [SingleR](https://github.com/dviraran/SingleR): (SingleR generate different cluster because it's annotated at single-cell level)\
    ![SingleR](data/processed/mice_integrated/mice_merged_SingleR_annotated.png)
