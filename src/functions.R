# This file contains the R functions that are used during the search for k-components

find_easy_2seps <- function(g) {
  # Find easy 2-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  numV <- length(V(g))

  exclude <- rep.int(0, numV)                   # Initialize filter
  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=3)               # Get all cliques of size >= 3
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    exclude[isolated] <- 1
  }

  nex <- length(exclude[which(exclude == 1)])

  cat("    Searching for easy 2-separators: #excluded: ", nex, "\n")

  for (e in E(g)) {
    # igraph version note
    # nodes <- ends(g,e)   # Use if version 1.0.0
    nodes <- get.edge(g,e) # Use if version 0.6.6

    if (exclude[nodes[1]] == 1 || exclude[nodes[2]] == 1) {
      next
    }
    result <- is.separator(g, nodes)
    if(result) {
      pc <- pc + 1
      pl[[pc]] <- sort(nodes)
    }
  }

  pl <- unique(pl)
  cat("    Found 2-separators: ", length(pl), "\n")
  return(pl)
}



find_all_2seps <-function(g) {
  # Find all 2-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=3)               # Get all cliques of size >= 3
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 2) {                 #   Keep track of 2-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for all 2-separators: #excluded: ", nex, "\n")

  for (i in 1:numV) {
    h <- gp - i
    arts <- articulation.points(h)
    for (art in sort(arts)) {
      pc <- pc + 1
      pl[[pc]] <- sort(c(V(gp)$orig[i], V(h)$orig[art]))
    }
  }

  pl <- unique(pl)
  cat("    Found 2-separators: ", length(pl), "\n")
  return(pl)
}



find_easy_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  numV <- length(V(g))

  exclude <- rep.int(0, numV)                   # Initialize filter
  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    exclude[isolated] <- 1
  }

  nex <- length(exclude[which(exclude == 1)])

  cat("    Searching for easy 3-separators: #excluded: ", nex, "\n")

  for (e in E(g)) {
    # igraph version note
    # nodes <- ends(g,e)   # Use if version 1.0.0
    nodes <- get.edge(g,e) # Use if version 0.6.6

    if (exclude[nodes[1]] == 1 || exclude[nodes[2]] == 1) {
      next
    }
    h <- g - c(nodes[1],nodes[2])
    arts <- articulation.points(h)
    for (art in sort(V(h)$orig[arts])) {
      pc <- pc + 1
      pl[[pc]] <- sort(c(nodes[1],nodes[2],art))
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}



find_all_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 3) {                 #   Keep track of 3-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for all 3-separators: #excluded: ", nex, "\n")

  for (i in 1:(numV-1)) {
    for (j in (i+1):numV) {
      nodes <- c(i,j)
      h <- gp  - nodes
      arts <- articulation.points(h)
      for (art in sort(arts)) {
        pc <- pc + 1
        pl[[pc]] <- sort(c(V(gp)$orig[i], V(gp)$orig[j], V(h)$orig[art]))
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}



find_tough_3seps <- function(g) {
  # Find all 3-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=4)               # Get all cliques of size >= 4
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 3) {                 #   Keep track of 3-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for tough 3-separators: #excluded: ", nex, "\n")

  psbl <- vector()
  for (v in V(gp)) {
    ns <- unlist(neighbors(gp, v))
    for (i in ns) {
      for (j in ns) {
        if(i > j && are.connected(gp,i,j) == 0) {
          vdp <- vertex.disjoint.paths(gp,i,j)
          if (vdp == 3) {
            psbl <- append(psbl, v)
          }
        }
      }
    }
  }
 
  psbl <- unique(psbl)
  cat("psbl: ", sort(psbl), "\n")
  np <- length(psbl)

  if (np > 0) {
    for (i in 1:(np-2)) {
      for (j in (i+1):(np-1)) {
        for (k in (j+1):np) {
          gpp <- gp - c(psbl[i], psbl[j], psbl[k])
          #result <- is.separator(gp, c(psbl[i], psbl[j], psbl[k]))
          if(!is.connected(gpp)) {
            pc <- pc + 1
            pl[[pc]] <- sort( c(V(gp)$orig[psbl[i]], V(gp)$orig[psbl[j]], V(gp)$orig[psbl[k]] ))
          }
        }
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 3-separators: ", length(pl), "\n")
  return(pl)
}



find_easy_4seps <- function(g) {
  # Find easy 4-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  numV <- length(V(g))

  exclude <- rep.int(0, numV)                   # Initialize filter
  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=5)               # Get all cliques of size >= 5
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    exclude[isolated] <- 1
  }

  nex <- length(exclude[which(exclude == 1)])
  cat("    Searching for easy 4-separators: #excluded: ", nex, "\n")

  nex <- length(exclude[which(exclude == 1)])
  gp <- g - exclude
  numE <- length(E(gp))

  for (i in 1:(numE-1)) {
    set1 <- get.edge(gp,E(gp)[i])
    cat("set1: ", set1, "\n")
    for (j in (i+1):numE) {
      set2 <- get.edge(gp,E(gp)[j])
      nodes <- c(set1[1], set1[2], set2[1], set2[2])
      nodes <- unique(nodes)
      if (length(nodes) == 4) {
        gpp <- gp - nodes
        if(!is.connected(gpp)) {
          pc <- pc + 1
          pl[[pc]] <- sort(c(V(gp)$orig[nodes[1]], V(gp)$orig[nodes[2]],
                             V(gp)$orig[nodes[3]], V(gp)$orig[nodes[4]]))
          cat(pl[[pc]], "\n")
        }
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 4-separators: ", length(pl), "\n")
  return(pl)
}



find_tough_4seps <- function(g) {
  # Find tough 4-vertex separators
  # Note - assumes that "orig" attribute has been set for graph

  pc   <- 0
  pl   <- list()
  excluded <- vector()

  d <- degree(g)                                # Recalculate node degrees
  cl <- maximal.cliques(g, min=5)               # Get all cliques of size >= 5
  for (i in 1:length(cl)) {                     # Loop over cliques
    clique <- unlist(cl[[i]])                   # Get members of clique
    dcl <- d[unlist(cl[[i]])]                   # Get degrees of clique members
    clsize <- length(clique)                    # Get size of clique
    worldly  <- clique[which(dcl >= clsize)]    # Get worldly nodes
    isolated <- clique[which(dcl <  clsize)]    # Get isolated nodes
    excluded <- append(excluded, isolated)      # Maintain list of included vertices
    if (length(worldly) == 4) {                 #   Keep track of 4-separators
      pc <- pc + 1                              #   since they may not be separators 
      pl[[pc]] <- sort(worldly)                 #   after excluding isolated nodes
    }
  }

  nex <- length(excluded)
  gp <- g - excluded
  numV <- length(V(gp))

  cat("    Searching for tough 4-separators: #excluded: ", nex, "\n")

  psbl <- vector()
  for (v in V(gp)) {
    #cat("    v: ", v, "\n") 
    ns <- unlist(neighbors(gp, v))
    for (i in ns) {
      for (j in ns) {
        if(i > j && are.connected(gp,i,j) == 0) {
          vdp <- vertex.disjoint.paths(gp,i,j)
          if (vdp == 4) {
            psbl <- append(psbl, v)
          }
        }
      }
    }
  }
 
  psbl <- unique(psbl)
  cat("psbl: ", sort(psbl), "\n")
  np <- length(psbl)

  if (np > 0) {
    for (i in 1:(np-3)) {
      for (j in (i+1):(np-2)) {
        for (k in (j+1):(np-1)) {
          for (l in (k+1):np) {
            gpp <- gp - c(psbl[i], psbl[j], psbl[k], psbl[l])
            #result <- is.separator(gp, c(psbl[i], psbl[j], psbl[k], psbl[l]))
            if(!is.connected(gpp)) {
              pc <- pc + 1
              pl[[pc]] <- sort( c(V(gp)$orig[psbl[i]], V(gp)$orig[psbl[j]], V(gp)$orig[psbl[k]], V(gp)$orig[psbl[l]] ))
	      cat(pl[[pc]], "\n")
            }
          }
        }
      }
    }
  }

  pl <- unique(pl)
  cat("    Found 4-separators: ", length(pl), "\n")
  return(pl)
}



kcore_delete_nodes <- function(g, kcomponent) {
  # Reduce a graph to k-cores by deleting nodes
  for (i in 1:100) {
    d <- degree(g)
    lowd_nodes <- V(g)[which(d < kcomponent)]
    g <- g - lowd_nodes
    if (length(lowd_nodes) == 0) {
      break
    }
  }
  return(g)
}



bicomponent_processing <- function(g, kcomponent, candidate) {
  # Find bicomponents
  bc <- biconnected.components(g)

  # Get size of largest bicomponent
  cutoff <- 0
  for(i in 1:bc$no) { 
    if(length(bc$components[[i]]) > cutoff) {
      cutoff <- length(bc$components[[i]])
    }
  }
  
  # Loop over bicomponents
  for (j in 1:bc$no) {
    gp <- induced.subgraph(g, bc$components[[j]])

    # Get size of bicomponent and skip if largest
    nverts_before <- length(V(gp))
    nedges_before <- length(E(gp))
    if(nverts_before == cutoff) {
      next
    }    

    ### Reduce to k-cores before writing out
    gp <- kcore_delete_nodes(gp, kcomponent)
    nverts <- length(V(gp))
    nedges <- length(E(gp))

    if(nverts <= kcomponent) {
    }
    if(nverts > kcomponent) {
      nedges_clique <- (nverts*(nverts-1))/2

      # Find size of minimum separator
      msep <- minimum.size.separators(gp)
      if(length(msep) > 0) {
        msep_min <- length(msep[[1]])
      } else {
        msep_min <- 0
      }

      if(nedges == nedges_clique) {
        candidate <- candidate + 1
	ksize <- length(V(gp))
        outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
        write.graph(gp, outfile, format="pajek")
      }
      else if(msep_min >= kcomponent) {
        candidate <- candidate + 1
        ksize <- length(V(gp))
        outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
        write.graph(gp, outfile, format="pajek")
      }
      else if(nedges != nedges_clique) {

        mwBlocks  <- cohesive.blocks(gp)
        blocks    <- blocks(mwBlocks)
        cohesion  <- cohesion(mwBlocks)
        parent    <- parent(mwBlocks)
        nblocks   <- length(blocks)

        for (ib in 1:nblocks) {
          if(cohesion[ib] >= kcomponent && cohesion[parent[ib]] < kcomponent) {
            gpp <- induced.subgraph(gp, blocks[[ib]])
            candidate <- candidate + 1
            ksize <- length(V(gpp))
            outfile <- paste("p", kcomponent, "_c", cohesion[ib], "_s", ksize, "_id", candidate, ".net", sep="")
            write.graph(gpp, outfile, format="pajek")
          }
        }

      }
    }
  }
  return(candidate)
}



reduce_to_largest_bicomponent <- function(g) {
  # Reduce graph to largest bicomponent
  bc <- biconnected.components(g)
  lbc <- vector(length=bc$no)
  for(i in 1:bc$no) { lbc[i] <- length(bc$components[[i]]) }
  maxbc <- which.max(lbc)
  g <- induced.subgraph(g, bc$components[[maxbc]])
  g <- set.vertex.attribute(g, "orig", V(g), V(g))
  return(g)
}



apply_separators <- function(g, candidate, sepsize, kcomponent, pl) {
  ### Use separators to reduce graph
  w <- g   # Make a copy of graph to work with
  for (i in 1:length(pl)) {
    nodes <- unlist(pl[[i]])
    h <- w - nodes
    clust <- clusters(h)
    ncl <- clust$no
    cutoff <- max(clust$csize) + sepsize

    for(j in 1:ncl) {
      setm <- (V(h)$orig[which(clust$membership == j)])
      setp <- union(nodes, setm)

      if(length(setp) == cutoff) {
        next
      }

      g[setm, V(g)] <- FALSE
      gp <- induced.subgraph(w, unlist(setp))

      ### Reduce to k-cores before writing out
      nverts_before <- length(V(gp))
      nedges_before <- length(E(gp))
      gp <- kcore_delete_nodes(gp, kcomponent)
      nverts <- length(V(gp))
      nedges <- length(E(gp))

      if(nverts <= kcomponent) {
      }
      if(nverts > kcomponent) {
        nedges_clique <- (nverts*(nverts-1))/2

        # Find size of minimum separator
        msep <- minimum.size.separators(gp)
        if(length(msep) > 0) {
          msep_min <- length(msep[[1]])
        } else {
          msep_min <- 0
        }

        if(nedges == nedges_clique) {
          candidate <- candidate + 1
          ksize <- length(V(gp))
          outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
          write.graph(gp, outfile, format="pajek")
        }
        else if(msep_min >= kcomponent) {
          candidate <- candidate + 1
          ksize <- length(V(gp))
          outfile <- paste("p", kcomponent, "_c", msep_min, "_s", ksize, "_id", candidate, ".net", sep="")
          write.graph(gp, outfile, format="pajek")
        }
        else if(nedges != nedges_clique) {
          mwBlocks  <- cohesive.blocks(gp)
          blocks    <- blocks(mwBlocks)
          cohesion  <- cohesion(mwBlocks)
          parent    <- parent(mwBlocks)
          nblocks   <- length(blocks)

          for (ib in 1:nblocks) {
            if(cohesion[ib] >= kcomponent && cohesion[parent[ib]] < kcomponent) {
              gpp <- induced.subgraph(gp, blocks[[ib]])
              candidate <- candidate + 1
              ksize <- length(V(gpp))
              outfile <- paste("p", kcomponent, "_c", cohesion[ib], "_s", ksize, "_id", candidate, ".net", sep="")
              write.graph(gpp, outfile, format="pajek")
            }
          }

        }
      }

    }
  }
  return(list(g, candidate))
}
