# Load necessary libraries
library(ggplot2)
library(pvclust)
library(ape)

# Read in the SSR data and class information
ssr = read.table("percentages_output.csv", sep = ",", header = TRUE, row.names = 1)
classes = read.table("Classes.txt", sep = "\t", header = TRUE)

# Filter out SSRs with zero standard deviation across samples
ssr_filter <- ssr[, apply(ssr, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# Sum of SSR content per sample
ssr_sum = data.frame(Sample = rownames(ssr_filter), SSR_Content = apply(ssr_filter, 1, sum))

# Principal Component Analysis (PCA) on filtered SSR data
pca_ssr_filter = prcomp(ssr_filter)
summary(pca_ssr_filter) # PC1 46.7%, PC2 11.8%

# Convert PCA results to data frame and merge with class information
pca_results = data.frame(Accession.ID = rownames(pca_ssr_filter$x), pca_ssr_filter$x)
pca_results = merge(pca_results, classes)

# PCA plot with labels
ggplot(data = pca_results, aes(x = PC1, y = PC2, col = Label)) + 
  geom_point(cex = 2) + 
  xlab("PC1 46.7%") + 
  ylab("PC2 11.8%")
ggsave("PCA_SSR_Labels.pdf", device = "pdf", width = 10, height = 10, dpi = 300)

# PCA plot with classes
ggplot(data = pca_results, aes(x = PC1, y = PC2, col = Class)) + 
  geom_point(cex = 2) + 
  xlab("PC1 46.7%") + 
  ylab("PC2 11.8%") + 
  theme(legend.position = "bottom")
ggsave("PCA_SSR_Classes.pdf", device = "pdf", width = 10, height = 10, dpi = 300)

# Write the distance matrix to an Excel file
write.table(x = as.matrix(dist(ssr_filter)), file = "Distance_SSR.xls", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Function to convert hierarchical clustering results into phylo format with node names
as.phylo.hclust.with.nodenames <- function (x, nodenames, ...) {
  # x: hierarchical clustering object
  # nodenames: Node names to include in the phylo object
  N <- dim(x$merge)[1]
  edge <- matrix(0L, 2 * N, 2)
  edge.length <- numeric(2 * N)
  node <- integer(N)
  node[N] <- N + 2L
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- x$merge[i, l]
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- x$height[i] - x$height[y]
      } else {
        edge[k, 2] <- -y
        edge.length[k] <- x$height[i]
      }
    }
    j <- j + 2L
  }
  if (is.null(x$labels)) x$labels <- as.character(1:(N + 1))
  node.lab <- nodenames[order(node)]
  obj <- list(edge = edge, edge.length = edge.length/2, tip.label = x$labels, Nnode = N, node.label = node.lab)
  class(obj) <- "phylo"
  reorder(obj)
}

# Read the distance matrix for SSR clustering
distance = read.table("distance_matrix_SSR_newLabels_all", sep = "\t", header = TRUE)
distance$New_Label = make.names(distance$Final_Label, unique = TRUE)
rownames(distance) = distance$New_Label
distance$New_Label = NULL
distance$Final_Label = NULL

# Perform hierarchical clustering with bootstrapping using pvclust
clustering = pvclust(parallel = TRUE, t(distance), method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = 1234)

# Extract bootstraps and convert clustering results to phylo format
bootstraps <- (round(clustering$edges, 2) * 100)[, 3]
yy <- as.phylo.hclust.with.nodenames(clustering$hclust, nodenames = bootstraps)

# Write the phylogenetic tree with bootstrap values to a Newick file
write.tree(phy = yy, file = "Bootstrap_HC_ALL.nwk")

# Plot tree with iTol (https://itol.embl.de/)