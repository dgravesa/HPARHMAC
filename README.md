HPAR HMAC

HPAR HMAC(High Performace Computing based Parallel Hiearchical Modal Association 
Clustering) is based on the Modalclust R package which performs Hierarchical Mode 
Association Clustering (HMAC) by http://www.scirp.org/journal/PaperInformation.aspx?PaperID=51471
Modal clustering techniques are especially designed to efficiently extract 
clusters in high dimensions with arbitrary density shapes. Further, clustering 
is performed over several resolutions and the results are summarized as a 
hierarchical tree, thus providing a model based multi resolution cluster analysis.

HPAR HMAC is an improvement over the parallel implementation in ModalCLust
as such it removes the limitation of the data size that can be processed as it 
supports MPI to enable the processing on multinode clusters and incorporates 
enhanced techniques to speed up the computation time.

