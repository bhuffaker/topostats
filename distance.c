/*
** Computes various distance-related graph metrics:
**
**   * average distance
**   * std deviation of distance
**
**   * average eccentricity
**   * graph radius
**   * graph diameter
**   * min degree in center
**   * min degree in periphery
**
** ---------------------------------------------------------------------
** Copyright (C) 2010 The Regents of the University of California.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <Judy.h>

#define false 0
#define true  1

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;
Pvoid_t index_node = (Pvoid_t)NULL;

Pvoid_t queue = (Pvoid_t)NULL;
Pvoid_t visited = (Pvoid_t)NULL;  /* Judy1 */
Pvoid_t distances = (Pvoid_t)NULL;
Pvoid_t distance_dist = (Pvoid_t)NULL;

Pvoid_t eccentricities = (Pvoid_t)NULL;
Pvoid_t distance_sum = (Pvoid_t)NULL;

unsigned long num_nodes, num_links;
unsigned long long num_pairs;  /* C(n, 2) pairs of nodes */

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
void dump_distance(char * filename);
void compute_distance_metrics(void);
unsigned long compute_node_distance_metrics(Word_t i0);
void compute_distance_statistics(void);


/* ====================================================================== */

void
load_graph(void)
{
  Word_t pi, i, v, l0, li, d, *pv;
  int rc;

  num_nodes = num_links = 0;

  pi = 0;  /* previous node ID; 0 == no previous node */
  l0 = 0;  /* starting index for the links of the current node */
  li = 1;  /* index of current link */
  d = 0;   /* node degree */

  while ((rc = scanf("%lu %lu", &i, &v)) == 2) {
    if (i == 0 || v == 0) {
      fputs("ERROR: node IDs must be > 0\n", stderr);
      exit(1);
    }

    if (i != pi) {
      if (pi != 0) {
	JLI(pv, links, l0);  *pv = d; /* save the degree of the previous node */
	l0 = li++;
	d = 0;
      }

      ++num_nodes;
      pi = i;
      JLI(pv, nodes, i);  *pv = l0;
      JLI(pv, index_node, num_nodes); *pv = i;
    }

    ++num_links;
    JLI(pv, links, li);  *pv = v;
    ++li;
    ++d;
  }

  if (rc != EOF && rc != 2) {
    fputs("ERROR: malformed line in input graph\n", stderr);
    exit(1);
  }

  if (num_nodes > 0) {
    JLI(pv, links, l0);  *pv = d;  /* save the degree of the last node */
  }

  num_pairs = (unsigned long long)num_nodes * (num_nodes - 1) / 2;
  num_links /= 2;  /* count undirected links */
  printf("loaded %lu nodes, %lu undirected links, %llu pairs\n",
	 num_nodes, num_links, num_pairs);
}


/* ====================================================================== */

void
dump_graph(void)
{
  Word_t i, *pv;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
    dump_links(i, *pv);
    JLN(pv, nodes, i);
  }
}


/* ====================================================================== */

void
dump_links(Word_t i, Word_t li)
{
  Word_t d, *pv;

  JLG(pv, links, li);
  d = *pv;

  while (d > 0) {
    --d;
    ++li;
    JLG(pv, links, li);
    printf("%lu %lu\n", i, *pv);
  }
}

/* ====================================================================== */

void 
dump_distance(char * filename) 
{
  Word_t i, li, deg, n, dist_sum, eccen, *pv;
  FILE * fp;


  fp = fopen(filename,"w");
  if (NULL == fp)
    fprintf (stderr, "failed to open file:%s\n", filename);

  printf("# storing distances in %s\n",filename);
  i = 1;
  fprintf(fp, "#node_id deg distance_avg eccentricity\n");
  JLG(pv, index_node, i); 
  while(NULL != pv) {
    n = *pv;
    JLG(pv, nodes, n);  
    if (NULL != pv) {
      li = *pv;
      JLG(pv, links, li);  deg = (NULL != pv) ? *pv : -1;
    }
    JLG(pv, distance_sum, n); dist_sum = (NULL != pv) ? *pv : -1;
    JLG(pv, eccentricities, n); eccen = (NULL != pv) ? *pv : -1;

    fprintf (fp, "%lu %lu %lf %lu\n", n, deg, 
	(1.0*dist_sum/(num_nodes-1)), eccen);
    JLG(pv, index_node, ++i);
  }
}

/* ====================================================================== */

/*
** Definitions:
**
**   eccentricity of node i = max distance from i
**   graph radius = min eccentricity in a graph
**   graph diameter = max eccentricity in a graph
**
**   graph periphery = nodes with max eccentricity
**   graph center = nodes with min eccentricity
*/
void
compute_distance_metrics(void)
{
  /* zero is never a valid distance or degree, so zero = sentinel */  
  unsigned long graph_radius=0, graph_diameter=0;
  unsigned long min_deg_in_graph_center=0, max_deg_in_graph_periphery=0;
  unsigned long eccentricity;
  double eccentricity_sum=0.0; /* for calculating avg eccentricity */
  Word_t i, li, deg, *pv;
  unsigned long num_nodes_in_graph_center = 0;
  unsigned long num_nodes_in_graph_periphery = 0;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
#ifdef DEBUG
    printf("\n* node %lu:\n", i);
#endif

    eccentricity = compute_node_distance_metrics(i);
    eccentricity_sum += eccentricity;

#ifdef DEBUG
    printf("  eccentricity %lu = %lu\n", i, eccentricity);
#endif

    /* find degree of the starting node */
    JLG(pv, nodes, i);  li = *pv;
    JLG(pv, links, li);  deg = *pv;

    if (eccentricity < graph_radius || graph_radius == 0) {
      printf("... decreasing radius from %lu to %lu with node %lu\n",
	     graph_radius, eccentricity, i);
      graph_radius = eccentricity;
      min_deg_in_graph_center = deg;
      num_nodes_in_graph_center = 1;
    }
    else if (eccentricity == graph_radius) {
      if ((deg < min_deg_in_graph_center || min_deg_in_graph_center == 0)) {
        min_deg_in_graph_center = deg;
      }
      num_nodes_in_graph_center++;
    }

    if (eccentricity > graph_diameter) {
      printf("... raising diameter from %lu to %lu with node %lu\n",
	     graph_diameter, eccentricity, i);
      graph_diameter = eccentricity;
      max_deg_in_graph_periphery = deg;
      num_nodes_in_graph_periphery = 1;
    }
    else if (eccentricity == graph_diameter) {
      if (deg > max_deg_in_graph_periphery) {
        max_deg_in_graph_periphery = deg;
      }
      num_nodes_in_graph_periphery++;
    }

    JLN(pv, nodes, i);
  }

  compute_distance_statistics();

  printf("average eccentricity = %.3f\n", eccentricity_sum / num_nodes);
  printf("graph radius = %lu\n", graph_radius);
  printf("graph diameter = %lu\n", graph_diameter);
  printf("min degree in graph center = %lu\n", min_deg_in_graph_center);
  printf("max degree in graph periphery = %lu\n", max_deg_in_graph_periphery);
  printf("num. nodes in graph center = %lu\n", num_nodes_in_graph_center);
  printf("num. nodes in graph periphery = %lu\n", num_nodes_in_graph_periphery);
}


/* ====================================================================== */

/*
** Computes and returns the eccentricity of the given node, and also updates
** the distance distribution.
**
** Distances are calculated with a breadth-first search.
*/
unsigned long
compute_node_distance_metrics(Word_t i0)
{
  Word_t n, qhead, qtail, i, i2, li, deg, dist, dist_sum, *pv, Rc_word;
  int Rc_int;
  unsigned long eccentricity=0;  /* zero is never a valid distance */

  dist_sum = 0;

  qhead = qtail = 0;
  JLI(pv, queue, qtail);  *pv = i0;  ++qtail;

  J1FA(Rc_word, visited);
  J1S(Rc_int, visited, i0);

  JLI(pv, distances, i0);  *pv = 0;

  while (qhead != qtail) {
    JLG(pv, queue, qhead);  i = *pv;
    JLD(Rc_int, queue, qhead);
    ++qhead;

    JLG(pv, distances, i);  dist = *pv;
    JLG(pv, nodes, i);  li = *pv;
    JLG(pv, links, li);  deg = *pv;

    while (deg > 0) {
      --deg;
      ++li;
      JLG(pv, links, li);  i2 = *pv;
#ifdef DEBUG
      printf("  ? %lu %lu\n", i, i2);
#endif

      J1T(Rc_int, visited, i2);
      if (!Rc_int) {
#ifdef DEBUG
	printf("  q %lu %lu = %lu\n", i, i2, dist + 1);
#endif
	J1S(Rc_int, visited, i2);
	JLI(pv, distances, i2);  *pv = dist + 1;
	JLI(pv, queue, qtail);  *pv = i2;  ++qtail;
	if (dist + 1 > eccentricity) {
	  eccentricity = dist + 1;
	}

	/* 
	** Used to store the average distance
	*/
	dist_sum += dist + 1;

	/*
	** Because we represent undirected links with a pair of directed links,
	** we need to ensure that we don't double count things.  We do this
	** by only gathering data when the link goes from a lower-numbered
	** node to a higher-numbered node.
	**
	**  NOTE: We must check against the *starting* node (i0) and not
	**        the current node being traversed (i).
	*/
	if (i0 < i2) {
#ifdef DEBUG
	  printf("  dist %lu %lu @%lu\n", i0, i2, dist + 1);
#endif
	  JLI(pv, distance_dist, dist + 1);  *pv += 1;
	}
      }
    }
  }

  JLI(pv, eccentricities, i0); *pv = eccentricity;
  JLI(pv, distance_sum, i0); *pv = dist_sum;

  return eccentricity;
}


/*
** Computes the average distance and the standard deviation of the
** distance distribution.
**
** The sample standard deviation code is efficient and numerically
** stable.  However, we compute the average distance separately using
** (mostly) integer arithmetic rather than simply re-using the mean
** calculated by the incremental standard deviation code because the
** latter seems to produce a mean that is slightly less accurate due
** to round-off errors.
*/
void
compute_distance_statistics(void)
{
  Word_t i, j, count, *pv;
  unsigned long long sum;  /* sum of the distances for avg distance calc */
  double sd_i, sd_mean, sd_q;

#ifdef DEBUG
  printf("\n* distance distribution:\n");
#endif

  sum = 0;
  sd_i = 0;
  sd_mean = sd_q = 0.0;

  i = 0;
  JLF(pv, distance_dist, i);
  while (pv != NULL) {
    count = *pv;
    sum += i * count;

#ifdef DEBUG
    printf("  %lu: %lu\n", i, count);
#endif

    /* std dev calculation */
    for (j = 0; j < count; j++) {
      sd_i += 1;
      sd_q += (sd_i - 1) * (i - sd_mean) * (i - sd_mean) / sd_i;
      sd_mean += (i - sd_mean) / sd_i;
    }

    JLN(pv, distance_dist, i);
  }

  printf("average distance = %.3f\n", sum / (double)num_pairs);
  printf("std deviation of distance = %.3f\n", sqrt(sd_q / (sd_i - 1)));
}


/* ====================================================================== */

int
main(int argc, char *argv[])
{
  int dump_distance_flag = false;
  int c;

  while ((c = getopt(argc, argv, "hd")) != -1) {
    switch (c) {
    case 'd':
      dump_distance_flag = true;
      break;

    case '?':
    case 'h':
    default:
      fprintf(stderr, "Usage: betweenness [-z]\n"
        "Options:\n"
	"  -h this message\n"
	"  -d dump to betweenness_link.txt and betweenness_nodes.txt\n"
        "  -z to normalize betweenness with n(n-1) [NOT recommended]\n");
      exit(1);
    }
  }
  load_graph();
#ifdef DEBUG
  dump_graph();
#endif
  compute_distance_metrics();

  if (true == dump_distance_flag) {
    dump_distance("distance_nodes.txt");
  }
  return 0;
}
