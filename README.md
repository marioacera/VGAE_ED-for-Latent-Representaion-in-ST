# VGAE_ED-for-Latent-Representaion-in-ST
Exploring Onthology based profiling for supporting raw mRNA counts information in Spatial Transcriptomics (ST) using a deep learning approach based on Variational Graph Autoencoders that follows a the guideline sets by SERD ( https://github.com/JinmiaoChenLab/SEDR ). This model poses an algorithm that is able to efficiently handle high dimensional data on clustering tasks, employing an iterative clustering loss as a soft spatial contrain. In this repository we pose ways to achieve alternative and refined latent representations of count data by exploiting deep learning and graph advantages.

Firstly, we present a complementary way of profiling cell or spots in Single Cell technologies using a Onthology reference, this information can be used alongside mRNA counts as support information or as base count information, depending on the sample and aims of the user.

To precisely capture tissue geometric borders we hypothesise that is important to account for non neighbourhood aggregated graph features. However when it comes to broadly capture anatomical relevant areas, using only the VGAE features as latent representation for clustering tasks can notably improve improve the performance of the method. Moreover, as it is shown below, Onthology based counts seems to better encode biological meaninfull features to capture tissue areas. 

Refined clustering results on the DLPFC dataset are shown below, just highlight that current state of the art models such as Bayes Space or SERD that achieve an ARI metric of 0.55 and 0.573.

 - Based on VGAE features from mRNA + Ontholgy counts: 0.578 ARI 

![show_Human_brain_leiden_regions_only_Graph_feat_noHVG](https://user-images.githubusercontent.com/56892292/134486112-9723f58a-7507-4f44-8753-7a2eb61040a9.png)![umap_Human_brain_leiden_regions_only_Graph_feat_noHVG](https://user-images.githubusercontent.com/56892292/136359556-b47925c0-662b-4241-bdc6-3dcb0b4c8e8d.png)


 - Based on VGAE features from Onthology counts: 0.629 ARI
![show_Human_brain_leiden_regions_only_Graph_feat_02_pc_dist_](https://user-images.githubusercontent.com/56892292/136359756-5fd1dc3c-9b7f-4a61-8b94-5d77d14106b0.png)
![umap_Human_brain_leiden_regions_only_Graph_feat_02_pc_dist_](https://user-images.githubusercontent.com/56892292/136359762-bc9dbee2-a608-4342-b929-15e5be876df6.png)

Seems reasonable to highlight that our model retrieve a latent representation, and from that point the clustering is performed following the state of the art pipeline in Single Cell (computing communities from the graph built based on the latent representation), therefore the UMAP is notably usefulll to asses the quality of the results.

On top of Onthology profiling, the flexibility that this model provides can be further exploited.
For example we can use a mixed metric to compute the adjacency matrix that define the graph structure in the VGAE. In order to help the model distinguishing tissue borders and transcriptomics gaps, we can bring in transcriptional information on the distance metric used to compute the input graph. It is shown below how the model accuracy improves when using a mixed distance metric jointly with Onthology based profiling:

 - Based on mixed metric VGAE features from Ontology counts: ARI 0.675
 
![show_Human_brain_leiden_regions_only_Graph_feat_02_pc_dist_](https://user-images.githubusercontent.com/56892292/136362335-db5a2808-6f5c-416e-a7a8-2cbe009aced1.png)
![umap_Human_brain_leiden_regions_only_Graph_feat_02_pc_dist_](https://user-images.githubusercontent.com/56892292/136362341-4e32097f-493e-4615-b293-160c57f9e2b1.png)




