require(igraph)

##### Parameters for k-components analysis
input_graph  <- "largest_quadconnected.net"
pre_output_graph <- "largest_prequintconnected.net"
output_graph <- "largest_quintconnected.net"
kcomponent   <- 5

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
  #cat("  Is connected: ", is.connected(g), "\n")
  if(lbc_size == lbc_size_old) break
  lbc_size_old <- lbc_size
}

cat("\nSTEP 2 - iteratively find/use 2-separators, reduce to bicomponents\n")
for (imaster in 1:100) {
  cat("---- STEP 2 iteration ", imaster, "\n")

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
  cat("  Apply 2-separators: Vertex count = ", size, "\n")

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  #cat("  Is connected: ", is.connected(g), "\n")
}

cat("\nSTEP 3 - Find/apply easy 3-separators\n")
for (imaster in 1:100) {
  cat("---- STEP 3 iteration ", imaster, "\n")

  ### Find and apply easy 3-separators
  pl <- find_easy_3seps(g)
  if(length(pl) == 0) break
  sepsize <- 3
  tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
  g         <- tmp[[1]]
  candidate <- tmp[[2]]

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply easy 3-separators: Vertex count = ", size, "\n")

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

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  #cat("  Is connected: ", is.connected(g), "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
    #cat("  Is connected: ", is.connected(g), "\n")
  }
}


# --------------------------
cat("\nSTEP 4 - Find/apply all 3-separators\n")
for (imaster in 1:100) {
  cat("---- STEP 4 iteration ", imaster, "\n")

  n3sep <- 0

  ### Find and apply tough 3-separators
  pl <- find_tough_3seps(g)
  n3sep <- n3sep + length(pl)
  if(length(pl) > 0) {
    sepsize <- 3
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply tough 3-separators: Vertex count = ", size, "\n")

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

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  #cat("  Is connected: ", is.connected(g), "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
    #cat("  Is connected: ", is.connected(g), "\n")
  }

  ### Find and apply easy 3-separators
  pl <- find_easy_3seps(g)
  n3sep <- n3sep + length(pl)
  if(length(pl) > 0) {
    sepsize <- 3
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }

  d <- degree(g)
  size <- length(V(g)[which(d > 0)])
  cat("  Apply easy 3-separators: Vertex count = ", size, "\n")

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

  # --- After k-coring, may have new lower degree separators ---

  ### Process bicomponents and save largest for further processing
  candidate <- bicomponent_processing(g, kcomponent, candidate)
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
  #cat("  Is connected: ", is.connected(g), "\n")

  ### Find and apply 2-separators
  for (ijk in 1:100) {
    pl <- find_all_2seps(g)
    if(length(pl) == 0) {break}
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]

    d <- degree(g)
    size <- length(V(g)[which(d > 0)])
    cat("  Apply 2-separators: Vertex count = ", size, "\n")
    
    ### Process bicomponents and save largest for further processing
    candidate <- bicomponent_processing(g, kcomponent, candidate)
    g <- reduce_to_largest_bicomponent(g)
    lbc_size <- length(V(g))
    cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")
    #cat("  Is connected: ", is.connected(g), "\n")
  }

  if(n3sep == 0) break

}

cat("\n\nSTEP 5 - Find/apply tough 4-separators\n")

for (ijk in 1:100) {
  ### Find and apply tough 4-separators
  pl <- find_tough_4seps(g)
  n4sep <- length(pl)
  if(length(pl) > 0) {
    sepsize <- 4
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  pl <- find_tough_3seps(g)
  n3sep <- length(pl)
  if(length(pl) > 0) {
    sepsize <- 3
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  pl <- find_all_2seps(g)
  n2sep <- length(pl)
  if(length(pl) > 0) {
    sepsize <- 2
    tmp       <- apply_separators(g, candidate, sepsize, kcomponent, pl)
    g         <- tmp[[1]]
    candidate <- tmp[[2]]
  }
  g <- reduce_to_largest_bicomponent(g)
  lbc_size <- length(V(g))
  cat("  Bicomponent processing: Vertex count = ", lbc_size, "\n")

  if(n4sep==0 && n3sep==0 && n2sep==0) {
    break
  }
}


write.graph(g, pre_output_graph, format="pajek")


cat("Now cleaning up with standard algorithm\n")
mydate <- date(); cat(mydate, "\n")


# Now cleanup the final graph
mwBlocks  <- cohesive.blocks(g)
blocks    <- blocks(mwBlocks)
cohesion  <- cohesion(mwBlocks)
parent    <- parent(mwBlocks)
nblocks   <- length(blocks)
for (ib in 1:nblocks) {
  if(cohesion[ib] >= kcomponent && cohesion[parent[ib]] < kcomponent) {
    gpp <- induced.subgraph(g, blocks[[ib]])
    candidate <- candidate + 1
    ksize <- length(V(gpp))
    outfile <- paste("p", kcomponent, "_c", cohesion[ib], "_s", ksize, "_id", candidate, ".net", sep="")
    write.graph(gpp, outfile, format="pajek")
  }
}


write.graph(g, output_graph, format="pajek")
mydate <- date(); cat(mydate, "\n")
