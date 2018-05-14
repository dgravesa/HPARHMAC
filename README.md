### High-Performance Accelerated Recursive Hierarchical Modal Association Clustering

*Posted to personal repository with permission of Dilip Patlolla and Sujith Nair*

*Publication may be found [here](https://www.osti.gov/biblio/1365649-high-performance-computing-based-parallel-hiearchical-modal-association-clustering-hpar-hmac)*

#### Overview

This package is based on the [Modalclust R](https://cran.r-project.org/web/packages/Modalclust/index.html) package, which is a multicore implementation of the [original HMAC](http://personal.psu.edu/jol2/hmac/) algorithm.
We reconsider the partitioning scheme used by the multicore R package and arrive at a solution that separates the MPI process layer from the sample partitioning layer.
The utility is developed with multicore and distributed systems in mind using MPI; however, the utility may also be built without dependence on MPI.

#### Improvements

By performing clustering with the minimum partition size parameter set, users will see improvement over the original runtime.
By dividing the data set into *N* partitions, the ideal speedup becomes *N*.
By using *M* processes in conjunction with partitioning, the ideal speedup becomes *MN*.
In our tests, the ideal speedup holds reasonably well for minimum partition sizes greater than 100.
For example, running a data set with two processes and a minimum partition size equal to one-tenth of the data size (thus 10 partitions) will result is approximately 20 times speedup.

#### Drawbacks

By using a minimum partition size parameter, there is inherent sampling in computing the density mixture.
This may result in some noise compared to the results of the original algorithm, where a minority of points may be misclustered.
Generally, these misclusterings occur on points with weak memberships -- they generally fall near cluster edges.
This same phenomenon occurs in the Modalclust R package.

Additionally, in our tests we found two types of clustering to occur.
In the first type, the number of clusters remains approximately constant as the number of points observed increases.
This is the ideal case for our algorithm, and is the typical case for clustering in general.
In the second type, the number of clusters increases approximately linearly as the number of points observed increases.
In our tests, this was most likely to occur with high-dimensional data, since significant difference in any dimension will cause points to be clustered separately in the HMAC algorithm.
Additionally, this may be causes by a bandwidth that is too low.
When the second type of result occurs, our algorithm performs more similarly to the original HMAC implementation.
This type of result is generally undesirable, but is a component of both the Mocalclust R package and original HMAC algorithm as well.
However, in cases where this *does* occur, our algorithm can still make use of multicore to process the top-level merge concurrently.

#### To Do

The Wiki is coming soon with some test case examples.
