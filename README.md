# VGAE_ED-for-Latent-Representaion-in-ST
Exploring Onthology based profiling for supporting raw mRNA counts information in Spatial Transcriptomics (ST) using a deep learning approach based on Variation Graph Autoencoders.

We present an alternative way of profiling cell or spots in Single Cell technologies using a Onthology reference, this information can be used alongside mRNA counts as support information or as base count information. 

Onthology profiling can help to capture anatomically meaninfull areas where mRNA count struggle to set the correct tissue borders.

Here are presented two results (VGAE_leiden plots, together with the ground truth annotated layers) that improve state of the art methods performance with an ARI metric of 0.578 and 0.588, respectively:


![show_Human_brain_leiden_regions_only_Graph_feat_noHVG](https://user-images.githubusercontent.com/56892292/134486112-9723f58a-7507-4f44-8753-7a2eb61040a9.png)
![show_Human_brain_leiden_regions_only_Graph_feat (copy)](https://user-images.githubusercontent.com/56892292/134486660-434a45d9-0eab-45aa-8a20-b5fdc12dfb89.png)

In contrast with two state of the art methods such as Bayes Space or SERD that achieve an ARI metric of 0.55 and 0.573.

