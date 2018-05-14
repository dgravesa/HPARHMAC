/*
	Daniel Graves
	HMAC.hpp
*/

#include <cmath>
#include <vector>
#include <map>
#include <list>

#define HMAC_MODE_THRES 1e-1
#define HMAC_DIST_THRES 1e-6

typedef std::vector<double> hmac_mode;

struct hmac_dataset
{
	std::vector<double>::const_iterator begin, end;
	unsigned long num_fields, num_coords;
};

/*
// return false if modes are sufficiently close
struct mode_lt
{
	bool operator()(const hmac_mode &left, const hmac_mode &right) const
	{
		for (int i = 0; i < left.size(); ++i)
		{
			if (left[i] < right[i] - HMAC_MODE_THRES) return true;
			if (left[i] > right[i] + HMAC_MODE_THRES) return false;
		}
		return false;
	}
};
*/

struct mode_lt
{
	bool operator()(const hmac_mode &left, const hmac_mode &right) const
	{
		for (int i = 0; i < left.size(); ++i)
			if (left[i] != right[i])
				return left[i] < right[i];

		return false;
	}
};

struct mode_eq
{
	bool operator()(const hmac_mode &left, const hmac_mode &right) const
	{
		for (int i = 0; i < left.size(); ++i)
			if (fabs(left[i] - right[i]) > HMAC_MODE_THRES)
				return false;

		return true;
	}
};

//typedef std::map<hmac_mode, std::list<unsigned long>, mode_lt> hmac_result; // TODO consider this heavily

// Recursive Modal Association Clustering
void RecursiveMAC(const hmac_dataset &dataset, double sigma, unsigned long min_part_size, std::vector<double> &modes, std::vector<unsigned long> &indices);
void RecursiveMAC(const hmac_dataset &dataset, double sigma, unsigned long min_part_size, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long offset = 0);

// Modal Association Clustering
void ModalAC(const hmac_dataset &dataset, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices);
void ModalAC(const hmac_dataset &dataset, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters);
void ModalAC(const hmac_dataset &dataset, const std::vector<double> &to_cluster, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices);
void ModalAC(const hmac_dataset &dataset, const std::vector<double> &to_cluster, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters);
void ModalAC(const hmac_dataset &dataset, const hmac_dataset &to_cluster, double sigma, std::vector<double> &modes, std::vector<unsigned long> &indices);
void ModalAC(const hmac_dataset &dataset, const hmac_dataset &to_cluster, double sigma, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long offset = 0);

// MAC output transforms
void Clusters2ModesAndIndices(const std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, unsigned long num_coords, std::vector<double> &modes, std::vector<unsigned long> &indices);
void ModesAndIndices2Clusters(const std::vector<double> &modes, const std::vector<unsigned long> &indices, unsigned long num_fields, std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters);

// Cluster map utilities
bool ClusterMapFind(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator &cluster);
void ClusterMapAdd(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, unsigned long index);
void ClusterMapMerge(std::map<hmac_mode, std::list<unsigned long>, mode_lt> &clusters, const hmac_mode &mode, std::list<unsigned long> &indices);

// Modal Expectation Maximization
void ModalEM(const hmac_dataset &dataset, double sigma, std::vector<double> &start, hmac_mode &mode);

