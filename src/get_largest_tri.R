require(igraph)

##### Parameters for k-components analysis
input_graph  <- "largest_biconnected.net"
output_graph <- "largest_triconnected.net"
kcomponent   <- 3

##### Include functions from "functions.R" file
source("functions.R")

mydate <- date(); cat(mydate, "\n")

g             <- read.graph(file=input_graph, format="pajek")
g             <- set.vertex.attribute(g, "orig", V(g), V(g))
candidate     <- 0 # Initialize k-candidate count
lbc_size_old  <- 0

cat("\nSTEP 1 - iteratively trim isolated nodes, reduce to k-cores, reduce to bicomponents\n")
for (imaster in 1:100) {
  cat("---- STEP 1 iteration ", imaster, "\n")

  ###### Process cliques with isolated nodes
  if(0) {
   # Skipping this step in search for tri-components. Generally results in longer run times
   for (iii in 1:100) {
     count <- 0                                    # Keep track of number of processed cliques
     d <- degree(g)                                # Recalculate node degrees
     cl <- maximal.cliques(g, min=kcomponent)      # Get all cliques of size >= kcomponent
     for (i in 1:length(cl)) {                     # Loop over cliques
       clique <- unlist(cl[[i]])                   # Get members of clique
       dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
       clsize <- length(clique)                    # Get size of clique
       isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
       worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes

       # If the clique contains exactly (k-1) worldly nodes, trim from main graph
       if(length(worldly) == (kcomponent-1) && length(isolated) > 0) {
         count <- count + 1

         # Write out the cliques that are not part of larger k-component
         if(clsize > kcomponent) { 
           # igraph version note
           #gp <- induced.subgraph(g, as_ids(clique)) # Use if version 1.0.0
           gp <- induced.subgraph(g, clique)          # Use if version 0.6.6
 	  msep_min <- length(V(gp)) - 1
           candidate <- candidate + 1
           ksize <- length(V(gp))
           outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
           write.graph(gp, outfile, format="pajek")
         }
         # igraph version note
         # g[as_ids(isolated), as_ids(clique)] <- FALSE # Use if version 1.0.0
         g[isolated, clique] <- FALSE                   # Use if version 0.6.6
       }

     }

     d <- degree(g)
     size <- length(V(g)[which(d > 0)])
     cat("  Clique processing: Vertex count = ", size, "\n")

     if(count == 0) break
   }
  }

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    # igraph version note
    #g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE # Use if version 1.0.0
    g[lowd_nodes, alld_nodes] <- FALSE                  # Use if version 0.6.6
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  if(lbc_size == lbc_size_old) break
  lbc_size_old <- lbc_size
}

cat("\nSTEP 2A - iteratively find/apply easy 2-separators, reduce to bicomponents\n")
for (imaster in 1:1) {
  cat("---- STEP 2A iteration ", imaster, "\n")

  ### Find 2-separators
  pl <- find_easy_2seps(g)
  if(length(pl) == 0) break

  ### Apply 2-separators
  sepsize <- 2
  tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
  g         <- tmp[[1]]
  candidate <- tmp[[2]]

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply easy 2-separators: Vertex count = ", size, "\n")

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    # igraph version note
    #g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE # Use if version 1.0.0
    g[lowd_nodes, alld_nodes] <- FALSE                  # Use if version 0.6.6
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
}

cat("\nSTEP 2B - iteratively find/apply all 2-separators, reduce to bicomponents\n")
for (imaster in 1:100) {
  cat("---- STEP 2B iteration ", imaster, "\n")

  ### Find 2-separators
  pl <- find_all_2seps(g)
  if(length(pl) == 0) break

  ### Apply 2-separators
  sepsize <- 2
  tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
  g         <- tmp[[1]]
  candidate <- tmp[[2]]

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply all 2-separators: Vertex count = ", size, "\n")

  ### Reduce to k-cores by deleting edges
  xold <- 0
  alld_nodes <- V(g)
  for (jjj in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    if (length(lowd_nodes) == xold) {
      break
    }
    xold <- length(lowd_nodes)
    # igraph version note
    #g[as_ids(lowd_nodes), as_ids(alld_nodes)] <- FALSE # Use if version 1.0.0
    g[lowd_nodes, alld_nodes] <- FALSE                  # Use if version 0.6.6
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  k-coring: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
}

write.graph(g, output_graph, format="pajek")
mydate <- date(); cat(mydate, "\n")
