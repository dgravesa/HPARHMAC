/*
	Daniel Graves
	HMAC.cpp
*/

#include "HMAC.hpp"

#include <cmath>
#include <algorithm>

#ifdef DEBUG
#include <iostream>
#endif

void RecursiveMAC(const hmac_dataset &dataset, double sigma, unsigned long min_part_size, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	std::map<hmac_mode, std::list<unsigned long>, mode_lt> clusters;

	// perform recursive MAC
	RecursiveMAC(dataset, sigma, min_part_size, clusters);

	// convert cluster map to modes and indices
	Clusters2ModesAndIndices(clusters, dataset.num_coords, modes, indices);
}

void RecursiveMAC(const hmac_dataset &dataset, double sigma, unsigned long min_part_size, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long offset)
{
	std::map<hmac_mode, std::list<unsigned long>, mode_lt> tmp_clusters;

	// default to base MAC
	if (min_part_size < 1 || min_part_size >= dataset.num_coords)
	{
		ModalAC(dataset, sigma, clusters);
		return;
	}

	// divide into left and right partitions and continue recursion
	if (dataset.num_coords > 2 * min_part_size)
	{
		// initialize partitions
		hmac_dataset left, right;
		left.num_coords = dataset.num_coords / 2;
		right.num_coords = dataset.num_coords - left.num_coords;
		left.num_fields = right.num_fields = dataset.num_fields;
		left.begin = dataset.begin;
		left.end = right.begin = dataset.begin + dataset.num_fields * left.num_coords;
		right.end = dataset.end;

		// execute left and right partitions
		RecursiveMAC(left, sigma, min_part_size, tmp_clusters, offset);
		RecursiveMAC(right, sigma, min_part_size, tmp_clusters, offset + left.num_coords);
	}

	// perform base MAC on windowed partitions
	else
	{
		// initialize to_cluster partitions
		hmac_dataset to_cluster_left, to_cluster_right;
		to_cluster_left.num_coords = dataset.num_coords / 2;
		to_cluster_right.num_coords = dataset.num_coords - to_cluster_left.num_coords;
		to_cluster_left.num_fields = to_cluster_right.num_fields = dataset.num_fields;
		to_cluster_left.begin = dataset.begin;
		to_cluster_left.end = to_cluster_right.begin = dataset.begin + dataset.num_fields * to_cluster_left.num_coords;
		to_cluster_right.end = dataset.end;

		// initialize windowed datasets
		hmac_dataset dataset_left, dataset_right;
		dataset_left.num_coords = dataset_right.num_coords = min_part_size;
		dataset_left.num_fields = dataset_right.num_fields = dataset.num_fields;
		dataset_left.begin = dataset.begin;
		dataset_left.end = dataset.begin + dataset.num_fields * dataset_left.num_coords;
		dataset_right.end = dataset.end;
		dataset_right.begin = dataset.end - dataset.num_fields * dataset_right.num_coords;

		// execute left and right partitions
		ModalAC(dataset_left, to_cluster_left, sigma, tmp_clusters, offset);
		ModalAC(dataset_right, to_cluster_right, sigma, tmp_clusters, offset + to_cluster_left.num_coords);
	}

	// perform merger by performing ModalEM on modes from tmp_clusters
	std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator tmp_cluster;
	for (tmp_cluster = tmp_clusters.begin(); tmp_cluster != tmp_clusters.end(); ++tmp_cluster)
	{
		std::vector<double> start(tmp_cluster->first), mode(dataset.num_fields);
		ModalEM(dataset, sigma, start, mode);
/*
		// get position to either newly inserted cluster or existing cluster
		std::list<unsigned long> dummy;
		std::pair<std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator, bool> modefind;
		modefind = clusters.insert(std::pair< hmac_mode, std::list<unsigned long> >(mode, dummy));

		// append membership from tmp_cluster to map cluster
		modefind.first->second.splice(modefind.first->second.end(), tmp_cluster->second);
*/
		ClusterMapMerge(clusters, mode, tmp_cluster->second);
	}
}

void ModalAC(const hmac_dataset &dataset, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	// perform MAC
	ModalAC(dataset, dataset, sigma, modes, indices);
}

void ModalAC(const hmac_dataset &dataset, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters)
{
	// perform MAC
	ModalAC(dataset, dataset, sigma, clusters);
}

void ModalAC(const hmac_dataset &dataset, const std::vector<double> &to_cluster, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	hmac_dataset to_cluster_dataset;
	to_cluster_dataset.begin = to_cluster.begin();
	to_cluster_dataset.end = to_cluster.end();
	to_cluster_dataset.num_fields = dataset.num_fields;
	to_cluster_dataset.num_coords = to_cluster.size() / dataset.num_fields;

	// perform MAC
	ModalAC(dataset, to_cluster_dataset, sigma, modes, indices);
}

void ModalAC(const hmac_dataset &dataset, const std::vector<double> &to_cluster, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters)
{
	hmac_dataset to_cluster_dataset;
	to_cluster_dataset.begin = to_cluster.begin();
	to_cluster_dataset.end = to_cluster.end();
	to_cluster_dataset.num_fields = dataset.num_fields;
	to_cluster_dataset.num_coords = to_cluster.size() / dataset.num_fields;

	// perform MAC
	ModalAC(dataset, to_cluster_dataset, sigma, clusters);
}

void ModalAC(const hmac_dataset &dataset, const hmac_dataset &to_cluster, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	std::map<hmac_mode, std::list<unsigned long>, mode_lt> clusters;

	// perform MAC
	ModalAC(dataset, to_cluster, sigma, clusters);

	// convert clusters map to modes and indices
	Clusters2ModesAndIndices(clusters, to_cluster.num_coords, modes, indices);
}

void ModalAC(const hmac_dataset &dataset, const hmac_dataset &to_cluster, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long offset)
{
	// build clustering with each point in to_cluster
	unsigned long to_cluster_index = offset;
	std::vector<double>::const_iterator coord;
	for (coord = to_cluster.begin; coord < to_cluster.end; coord += dataset.num_fields)
	{
		std::vector<double> start(coord, coord + dataset.num_fields), mode(dataset.num_fields);
		ModalEM(dataset, sigma, start, mode);
/*
		std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator modefind = clusters.find(mode);

		// existing mode found
		if (modefind != clusters.end()) modefind->second.push_back(to_cluster_index++); // TODO test merge used in recursive MAC using std c++11

		// insert new mode
		else
		{
			std::list<unsigned long> dummy{to_cluster_index++};
			clusters.insert(std::pair< hmac_mode, std::list<unsigned long> >(mode, dummy));
		}
*/
		ClusterMapAdd(clusters, mode, to_cluster_index++);
	}
}

void Clusters2ModesAndIndices(const std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long num_coords, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	unsigned long num_fields = clusters.begin()->first.size();

	// allocate output vectors
	modes.resize(num_fields * clusters.size());
	indices.resize(num_coords);

	unsigned long mode_index = 0;
	std::map<hmac_mode, std::list<unsigned long>, mode_lt>::const_iterator cluster;
	for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster, ++mode_index)
	{
		// copy mode into vector
		std::copy_n(cluster->first.begin(), num_fields, modes.begin() + num_fields * mode_index);

		// set coordinate mapping to mode index
		std::list<unsigned long>::const_iterator coord_index;
		for (coord_index = cluster->second.begin(); coord_index != cluster->second.end(); ++coord_index)
			indices[*coord_index] = mode_index;
	}
}

void ModesAndIndices2Clusters(const std::vector<double> &modes, const std::vector<unsigned long> &indices, unsigned long num_fields, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters)
{
	unsigned long num_coords = indices.size();

	// initialize map with modes
	std::vector<double>::const_iterator mode;
	for (mode = modes.begin(); mode < modes.end(); mode += num_fields)
	{
		hmac_mode m(mode, mode + num_fields);
		std::list<unsigned long> dummy;
		clusters.insert(std::pair< hmac_mode, std::list<unsigned long> >(m, dummy));
	}

	// create vector of cluster locations in map
	unsigned long mode_index = 0;
	std::vector<std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator> map_clusters(clusters.size());
	std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator cluster;
	for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster)
		map_clusters[mode_index++] = cluster;

	// add coordinate indices to mode membership lists
	for (unsigned long coord_index = 0; coord_index < num_coords; ++coord_index)
		map_clusters[indices[coord_index]]->second.push_back(coord_index);
}

bool ClusterMapFind(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator &cluster)
{
	hmac_mode low(mode), high(mode);

	// initialize low and high values
	for (int i = 0; i < mode.size(); ++i)
	{
		low[i] -= HMAC_MODE_THRES;
		high[i] += HMAC_MODE_THRES;
	}

	// get lower and upper bounds
	std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator low_bound, high_bound;
	low_bound = clusters.lower_bound(low);
	high_bound = clusters.upper_bound(high);

	// search bounds for equivalent mode key
	mode_eq meq;
	for (cluster = low_bound; cluster != high_bound; ++cluster)
		if (meq(cluster->first, mode))
			return true;

	return false;
}

void ClusterMapAdd(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, unsigned long index)
{
	std::list<unsigned long> dummy{index};
	ClusterMapMerge(clusters, mode, dummy);
}

void ClusterMapMerge(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, std::list<unsigned long> &indices)
{
	// search for cluster in map
	bool found;
	std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator cluster;
	found = ClusterMapFind(clusters, mode, cluster);

	// cluster not found, insert new cluster
	if (!found)
	{
		std::list<unsigned long> dummy;
		cluster = clusters.insert(cluster, std::pair< hmac_mode, std::list<unsigned long> >(mode, dummy));
	}

	// add indices to cluster
	cluster->second.splice(cluster->second.end(), indices);
}

void ModalEM(const hmac_dataset &dataset, double sigma, std::vector<double> &start, hmac_mode &mode)
{
	//unsigned long num_coords = (data_end - data_begin) / num_fields;
	std::vector<double> p(dataset.num_coords);
	double sigma_sq = sigma * sigma;
	double dist_sq;

	do
	{
		double p_sum = 0.0;

		// calculate density impacts from each coordinate
		for (unsigned long i = 0; i < dataset.num_coords; ++i)
		{
			std::vector<double>::const_iterator coord = dataset.begin + dataset.num_fields * i;

			dist_sq = 0.0;
			for (unsigned long j = 0; j < dataset.num_fields; ++j)
				dist_sq += (coord[j] - start[j]) * (coord[j] - start[j]);

			p[i] = exp(-0.5 * dist_sq / sigma_sq);
			p_sum += p[i];
		}

		// normalize p
		for (unsigned long i = 0; i < dataset.num_coords; ++i)
			p[i] /= p_sum;

		for (unsigned long j = 0; j < dataset.num_fields; ++j)
			mode[j] = 0.0;

		// calculate new mode position
		for (unsigned long i = 0; i < dataset.num_coords; ++i)
		{
			std::vector<double>::const_iterator coord = dataset.begin + dataset.num_fields * i;

			for (unsigned long j = 0; j < dataset.num_fields; ++j)
				mode[j] += p[i] * coord[j];
		}

		dist_sq = 0.0;
		for (unsigned long j = 0; j < dataset.num_fields; ++j)
			dist_sq += (mode[j] - start[j]) * (mode[j] - start[j]);

		start = mode;
	}
	while (dist_sq > HMAC_DIST_THRES * sigma);
}

