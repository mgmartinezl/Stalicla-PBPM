# Title     : Clustering PTV patients
# Objective : Extract clusters of PTV patients with hierarchical techniques
# Created by: Gabriela Martinez - airamgabriela17@gmail.com
# Created on: 9/6/19

# The input dataset for this exercise is the binary matrix generated
# from the PBPM protocol. In this particular run, only PTV patients were considered.

#----- Data reading and processing ------#

path = "/home/gabi/Documents/CookieCutter/staliclaProjectTemplate/PBPM/data/raw/binary-matrix-2019-08-23_13꞉08꞉35.csv"

df = read.csv(path, header=TRUE, sep=",", skip = 1)

# Drop features with only 0s or 1s
df <- df[, colSums(df != 0) > 0]

#----- Similarity matrix for patients -----#

# Compute similarity matrix for patients, using Jaccard's distance for binary dichotomic variables
patients <- df[,2:ncol(df)]
patients.simmilarity <- dist(patients, method="binary")
patients.simmilarity.matrix <- as.matrix(patients.simmilarity)

#----- Hierarchical clustering for patients -----#

# Compute clusters
library(FactoClass)
patients.clustering <- ward.cluster(patients.simmilarity, plots=TRUE, h.clust = 1)

# Visualize dendrogram
library(ggdendro)
ggdendrogram(patients.clustering)

# Extract clusters >> 2 to 5
patients.clusters <- cutree(patients.clustering, k=2:5)

# Append clusters to patient IDs
library(dplyr)
patients.clusters <- mutate(df,
                           C2 = patients.clusters[,1],
                           C3 = patients.clusters[,2],
                           C4 = patients.clusters[,3],
                           C5 = patients.clusters[,4])

#Count the occurrences for three clusters
#count(patients.clusters, C3)

# Activate this to only extract IDs and clusters (not pathways)
patients.clusters <- dplyr::select(patients.clusters, child_id, C2, C3, C4, C5)
write.csv(patients.clusters, file = "analysis/PTV-patients-hc.csv", row.names=FALSE)

#----- Similarity matrix for pathways -----#

# Compute phi coefficient of correlation for pathways
pathways.correlation <- cor(patients, use = "pairwise.complete.obs")

# Compute similarity matrix for pathways
pathways.simmilariy = 1 - pathways.correlation
pathways.distance = as.dist(pathways.simmilariy)

#----- Hierarchical clustering for pathways -----#

# Compute clusters
library(FactoClass)
pathways.clustering <- ward.cluster(pathways.distance, plots=TRUE, h.clust = 1)

# Visualize dendrogram
library(ggdendro)
ggdendrogram(pathways.clustering)

# Extract clusters >> 2 to 5
pathways.clusters <- cutree(pathways.clustering, k=2:5)

# Activate this to only extract pathways and clusters (not pathways)
write.csv(pathways.clusters, file = "analysis/PTV-pathways-hc.csv", row.names=TRUE)
