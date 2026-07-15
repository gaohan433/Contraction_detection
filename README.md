\# Contraction Detection



This repository supports the manuscript:



"Enhancing Uterine Contraction Detection in Electromyometrial Imaging Based on Spatial Continuity of the Electric Field"
> https://doi.org/10.1016/j.bspc.2026.110537


The repository contains the contraction-detection methods evaluated in the study, MATLAB scripts for statistical analysis and plotting, data, and generated results.



\---



\## 📁 Repository Structure



\### `Methods\_except\_for\_ts2vec/`



This folder contains the implementations of all contraction-detection methods evaluated in the study, except for TS2Vec.



The sample entropy method is based on the contraction-detection approach described in:



> Chen Z, Wang M, Zhang M, Huang W, Feng Y, Xu J.  

> "Automatic detection and characterization of uterine contraction using electrohysterography."  

> Biomedical Signal Processing and Control. 2024.  

> https://doi.org/10.1016/j.bspc.2023.105840



The implementation included in this repository was adapted from the authors’ publicly available code:



https://github.com/zhenqinchen/EHG-contraction



Modifications were made as needed to apply the method to the EMMI/EHG data and evaluation framework used in the accompanying study. Please refer to the original paper and repository for complete methodological details and applicable licensing terms.



\### `ts2vec/`



This folder contains the TS2Vec-based contraction-detection pipeline.



The implementation was developed based on the original TS2Vec repository:



https://github.com/zhihanyue/ts2vec



TS2Vec was introduced in:



> Yue Z, Wang Y, Duan J, Yang T, Huang C, Tong Y, Xu B.  

> "TS2Vec: Towards Universal Representation of Time Series."  

> Proceedings of the AAAI Conference on Artificial Intelligence. 2022.

> https://arxiv.org/abs/2106.10466

Please refer to the original repository for additional implementation details and applicable licensing terms.



\### `Matlab/`



This folder contains the MATLAB scripts used for:



\- Statistical analysis

\- Performance comparison

\- Result visualization

MATLAB was used for statistical analysis and plotting in this study.



\### `data/`



This folder contains data required to run the available analyses.



Any data included in this repository are provided only for research and reproducibility purposes. Restricted, identifiable, or protected participant data are not distributed through this repository.



\### `result/`



This folder contains generated outputs, intermediate results, statistical summaries, and files used for visualization.



\---



\## 🧠 Code Overview



The contraction-detection methods are organized into two main groups:



1\. Conventional and feature-based methods implemented in `Methods\_except\_for\_ts2vec/`

2\. The TS2Vec-based method implemented in `ts2vec/`



The outputs from these methods are evaluated statistically and visualized using the scripts provided in `Matlab/`.



