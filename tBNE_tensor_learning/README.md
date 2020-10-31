# Brain Network Embedding

The task of converting brain network data from graph structures to vectorial representations is formulated as a brain network embedding problem. We leverage tensor factorization techniques to obtain latent representations of brain networks which can be further utilized to facilitate downstream tasks. In particular, undirected brain networks are stacked as a partially symmetric tensor before conducting factorization. The self-report data are incorporated as guidance in the tensor factorization procedure to learn latent factors that are consistent with the side information. Furthermore, the representation learning and classifier training are blended into a unified optimization framework to obtain discriminative representations, by allowing the classifier parameters to interact with the original brain network data via latent factors and the representation learning process to be aware of the supervision information. The formulated optimization problem can be interpreted as partially coupled matrix and tensor factorization with constraints.

License
-------
© Bokai Cao, 2018. Licensed under an [Apache-2](https://github.com/caobokai/brain-network-embedding/blob/master/LICENSE) license.

Reference
---------
Bokai Cao, Lifang He, Xiaokai Wei, Mengqi Xing, Philip S. Yu, Heide Klumpp and Alex D. Leow. [t-BNE: Tensor-based Brain Network Embedding](https://www.cs.uic.edu/~bcao1/doc/sdm17.pdf). In SDM 2017.
