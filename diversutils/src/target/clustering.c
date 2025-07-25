#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <stdint.h>

#include "graph.h"
#include "distances.h"

int32_t davies_bouldin_index(struct graph* const g, const int32_t* const cluster_index_per_node, float* res){
	int32_t num_clusters;
	uint64_t i, j;
	float result, local_maximum, r_ij;
	float* centroids, *cluster_diameters;
	int32_t* num_elements_per_cluster;
	size_t alloc_size;

	num_clusters = cluster_index_per_node[0];

	// find number of clusters
	for(i = 1 ; i < g->num_nodes ; i++){
		if(cluster_index_per_node[i] > num_clusters){
			num_clusters = cluster_index_per_node[i];
		}
	}

	num_clusters++; // because of base 0

	if(num_clusters <= 0){
		perror("num_clusters <= 0\n");
		return 1;
	}

	// printf("num_clusters: %i\n", num_clusters);

	alloc_size = num_clusters * g->num_dimensions * sizeof(float);
	centroids = malloc(alloc_size);
	if(centroids == NULL){
		perror("failed to malloc\n");
		return 1;
	}
	memset(centroids, '\0', alloc_size);

	alloc_size = num_clusters * sizeof(int32_t);
	num_elements_per_cluster = malloc(alloc_size);
	if(num_elements_per_cluster == NULL){
		perror("failed to malloc\n");
		free(centroids);
		return 1;
	}
	memset(num_elements_per_cluster, '\0', alloc_size);

	alloc_size = num_clusters * sizeof(float);
	cluster_diameters = malloc(alloc_size);
	if(num_elements_per_cluster == NULL){
		perror("failed to malloc\n");
		free(centroids);
		free(num_elements_per_cluster);
		return 1;
	}
	memset(cluster_diameters, '\0', alloc_size);

	// compute centroids
	for(i = 0 ; i < g->num_nodes ; i++){
		if(cluster_index_per_node[i] < 0){continue;}
		for(j = 0 ; j < (uint64_t) g->num_dimensions ; j++){
			centroids[cluster_index_per_node[i] * g->num_dimensions + j] += g->nodes[i].vector.fp32[j];
		}
		num_elements_per_cluster[cluster_index_per_node[i]] += 1;
	}

	for(i = 0 ; i < (uint64_t) num_clusters ; i++){
		for(j = 0 ; j < (uint64_t) g->num_dimensions ; j++){
			centroids[i * g->num_dimensions + j] /= (float) num_elements_per_cluster[i];
		}
	}

	// compute cluster diameters
	for(i = 0 ; i < g->num_nodes ; i++){
		if(cluster_index_per_node[i] < 0){continue;}
		/*
		for(j = 0 ; j < g->num_nodes ; j++){
			if(cluster_index_per_node[i] != cluster_index_per_node[j]){continue;}
			cluster_diameters[cluster_index_per_node[i]] += cosine_distance_fp32(g->nodes[i].vector.fp32, g->nodes[j].vector.fp32, g->num_dimensions);
		}
		*/
		cluster_diameters[cluster_index_per_node[i]] += cosine_distance_fp32(g->nodes[i].vector.fp32, &(centroids[cluster_index_per_node[i] * g->num_dimensions]), g->num_dimensions);
	}
	for(i = 0 ; i < (uint64_t) num_clusters ; i++){
		cluster_diameters[i] /= (float) num_elements_per_cluster[i];
	}

	// compute DB-index
	result = 0.0f;
	for(i = 0 ; i < (uint64_t) num_clusters ; i++){
		// local_maximum = -1.0f;
		local_maximum = 0.0f;
		for(j = 0 ; j < (uint64_t) num_clusters ; j++){
			if(i == j){continue;}
			r_ij = (cluster_diameters[i] + cluster_diameters[j]) / cosine_distance_fp32(&(centroids[i * g->num_dimensions]), &(centroids[j * g->num_dimensions]), g->num_dimensions);
			if(r_ij > local_maximum){
				local_maximum = r_ij;
			}
			// printf("r_ij %lu %lu %f\n", i, j, r_ij);
		}
		result += local_maximum;
	}
	result /= (float) num_clusters;

	(*res) = result;

	free(cluster_diameters);
	free(num_elements_per_cluster);
	free(centroids);

	return 0;
}

#endif
