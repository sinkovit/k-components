#H2 k-components identification software

This directory contains the source code for generating the largest 1-,
2-, 3-, 4- and 5-components of a graph.

* Code requires the igraph library. Currently compatible with igraph
  0.6.6, but commented code included to allow use of igraph
  1.0.0. Search for comments containing the line "igraph version note"

* Currently assumes Pajek format, but this could proabably be changed
  easily since igraph can also handle other graph formats

* This is currently research code, rather than production quality
  software. We will be improving the code over time

* The software may have some minor issues as we reach the final
  stages. For example, in the processing of the "social.net" data set,
  the largest 5-component is actually trimmed off of the main
  graph. This may require manual examination of the p_*.net files
  (described below)

Source code files are described below

#H3 functions.R

Collection of R functions that are used by multiple scripts

#H3 get_largest_bi.R

For convenience, script takes no command line arguments but does
require that the starting graph be named "start.net". Will produce two
new files named "largest_connected.net" and "largest_biconnected.net".

#H3 get_largest_tri.R

Generates "largest_triconnected.net" starting from "largest_biconnected.net".

#H3 get_largest_quad.R

Generates "largest_quadconnected.net" starting from "largest_triconnected.net".

#H3 get_largest_quint.R

Generates "largest_quintconnected.net" starting from
"largest_quadconnected.net". Also writes out a file named
"largest_prequintconnected.net" that can be used as a starting point
for completing search for higher-order k-components using the standard
max-flow algorithms.

#H3 unique_candidates.pl

Script used to remove duplicates of the small components generated
during the search for the largest k-components. To use, execute the
following sequence of commands in the directory containg the results

* mkdir UNIQUE

* md5sum p*_c*_s*.net > digests

* ./unique_candidates.pl digests

#H2 Additional notes

The "get_largest_*.R" scripts also generate smaller components that
were trimmed from the primary graph during the graph reduction phases
and application of vertex separators. These graphs are named using
the convention

p[A]_c[B]_s[C]_id[X].net

where,

* A = the k-component that is being searched for

* B = connectivity of the component

* C = size (number of vertices) of component

* X = numerical identifier for the component (used primarily to
  distinguish components with the same values for A, B, and C)
