# VGAE_ED-for-Latent-Representaion-in-ST
Exploring Onthology based profiling for supporting raw mRNA counts information in Spatial Transcriptomics (ST) using a deep learning approach based on Variational Graph Autoencoders that follows a the guideline sets by SERD ( https://github.com/JinmiaoChenLab/SEDR ). This model poses an algorithm that is able to efficiently handle high dimensional data on clustering tasks, employing an iterative clustering loss as a soft spatial contrain. 

We present an alternative way of profiling cell or spots in Single Cell technologies using a Onthology reference, this information can be used alongside mRNA counts as support information or as base count information, depending on the sample and aims of the user 

We prove here that Onthology profiling can help to capture anatomically meaninfull areas where mRNA count struggle to set the correct tissue borders.

Here are presented two results (VGAE_leiden plots, together with the ground truth annotated layers) that improve state of the art methods performance with an ARI metric of 0.578(using mRNA and Ontology counts) and 0.606, respectively:


![show_Human_brain_leiden_regions_only_Graph_feat_noHVG](https://user-images.githubusercontent.com/56892292/134486112-9723f58a-7507-4f44-8753-7a2eb61040a9.png)
![show_Human_brain_leiden_regions_only_Graph_feat_02_pc_dist_](https://user-images.githubusercontent.com/56892292/136046622-21c32247-3641-4288-abec-a829716d8bf5.png)

In contrast with two state of the art methods such as Bayes Space or SERD that achieve an ARI metric of 0.55 and 0.573.

