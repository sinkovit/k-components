require(igraph)

# Read graph and get rid of multiple edges between nodes
g <- read.graph("start.net", format="pajek")
g <- simplify(g)

# Reduce to largest connected cluster and save
cat("Number of clusters: ", clusters(g)$no, "\n\n")
cat("-- Cluster sizes --\n")
# print(sort(clusters(g)$csize))
maxcl <- which.max(clusters(g)$csize)                     # Identify largest cluster
subv <- V(g)[ which(clusters(g)$membership == maxcl) ]    # Get members of largest cluster
g <- induced.subgraph(g, subv)                            # Replace g with largest cluster
write.graph(g, "largest_connected.net", format="pajek")

# Extract largest biconnected component and save
bc <- biconnected.components(g)
lbc <- vector(length=length(bc$components))
for(i in 1:length(bc$components)) { lbc[i] <- length(bc$components[[i]]) }
cat("Number of biconnected components: ", bc$no, "\n\n")
cat("-- biconnected component sizes --\n")
# print(sort(lbc))
maxbc <- which.max(lbc)
g <- induced.subgraph(g, bc$components[[maxbc]])
write.graph(g, "largest_biconnected.net", format="pajek")
