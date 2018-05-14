/*
	Daniel Graves
	HMACTool.cpp
*/

#include "HMAC.hpp"

#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef WITH_BOOST
#include <boost/filesystem.hpp>
#endif

//#include <getopt.h>
#include "Opts.hpp"

#include <cstdlib>

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <vector>
#include <random>
#include <algorithm>
#include <utility>
#include <chrono>

/*
static struct option long_options[] =
{
	{ "help", no_argument, NULL, 'h' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "quiet", no_argument, NULL, 'q' },
	{ "sigmas", required_argument, NULL, 's' },
	{ "recursive", required_argument, NULL, 'r' },
	{ "normalize", optional_argument, NULL, 'N' },
	{ "shuffle", optional_argument, NULL, 'S' },
	{ "resume", required_argument, NULL, 'R' },
	{ "output-dir", required_argument, NULL, 'o' },
	{ "output-format", required_argument, NULL, 'f' },
	{ 0, 0, 0, 0 }
};
*/

static Option options[] =
{
	{ "help", 'h', NO_ARG, "print this usage dialog" },
	{ "verbose", 'v', NO_ARG, "print verbose runtime information to stdout" },
	{ "quiet", 'q', NO_ARG, "suppress output to stdout" },
	{ "sigmas", 's', REQUIRED_ARG, "specify sigmas for computation" },
	{ "recursive", 'r', REQUIRED_ARG, "perform recursive HMAC using minimum partition size specified" },
	{ "normalize", 'N', OPTIONAL_ARG, "use data normalization for computation (default VALUE=yes)" },
	{ "shuffle", 'S', OPTIONAL_ARG, "use data shuffling for computation (default VALUE=yes)" },
	{ "resume", 'R', REQUIRED_ARG, "specify most recent vectorized output to continue hierarchy construction" },
	{ "output-dir", 'o', REQUIRED_ARG, "set directory to write outputs for each sigma (default VALUE=./)" },
	//{ "output-format", 'f', REQUIRED_ARG, "set format of output files [(vectorized), clusters, none]" },
	{ 0, 0, 0, 0 }
};

////static const char short_options[] = "hv::s:r:N::S::o:n:f:";
//static const char short_options[] = "hvqs:r:N::S::R:o:f:";

bool open_dataset(const std::string &file_name, std::vector<double> &coords, unsigned long &num_coords, unsigned long &num_fields);
bool open_vectorized(const std::string &file_name, std::vector<double> &modes, std::vector<unsigned long> &indices, unsigned long &num_modes, unsigned long &num_coords, unsigned long &num_fields);

#ifdef WITH_MPI
void SendModesAndIndices(int to_proc_id, std::vector<double> &modes, std::vector<unsigned long> &indices);
void RecvModesAndIndices(int from_proc_id, std::vector<double> &modes, std::vector<unsigned long> &indices);

void hmac_exit(int code)
{
	MPI_Finalize();
	exit(code);
}
#else
void hmac_exit(int code)
{
	exit(code);
}
#endif

// default, only changes when using MPI
int num_procs = 1, proc_id = 0;

int main(int argc, char **argv)
{
	std::string input_fname;
	std::list<double> sigmas;
	unsigned long min_part_size = 0;
	bool verbose = false, quiet = false, normalize = true, shuffle = true, resume = false;
	std::string output_dir = ".", output_name;
	std::string resume_input_fname;
	enum hmac_output_format { hof_clusters, hof_vectorized, hof_none } output_format = hof_clusters;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
#endif

	// parse command line options
	OptHandler opt_handler(options, argc, argv);
	//int opt_code;
	char opt_code;
	std::string opt_arg;
	//while ((opt_code = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
	while (opt_handler.getOpt(opt_code, opt_arg))
	{
		//if (optarg != NULL) opt_arg = optarg;
		//else opt_arg.clear();
		std::stringstream ss(opt_arg);

		switch (opt_code)
		{
			// set verbose flag
			case 'v':
				verbose = true;
				break;

			// set quiet flag
			case 'q':
				quiet = true;
				break;

			// set sigmas for calculation
			case 's':
				while (!ss.eof())
				{
					char c = ss.peek();
					if (isdigit(c) || c == '-' || c == '.')
					{
						// get sigma values
						double val;
						ss >> val;

						// test value
						if (val <= 0.0)
						{
							if (proc_id == 0) std::cerr << "error: sigmas must be greater than 0.0\n";
							hmac_exit(5);
						}
						else if (sigmas.size() > 0 && val <= sigmas.back())
						{
							if (proc_id == 0) std::cerr << "error: sigmas must be an increasing list\n";
							hmac_exit(6);
						}

						// add sigma to list
						sigmas.push_back(val);
					}
					else ss.get();
			}
				break;

			// set recursion parameter
			case 'r':
				//min_part_size = std::stoul(optarg);
				min_part_size = std::stoul(opt_arg.c_str());
				break;

			// set or unset normalization flag
			case 'N':
				if (opt_arg.empty() || opt_arg == "yes") normalize = true;
				else if (opt_arg == "no") normalize = false;
				else if (proc_id == 0) std::cerr << "warning: normalize option expects 'yes' or 'no'\n";
				break;

			// set or unset data shuffle flag
			case 'S':
				if (opt_arg.empty() || opt_arg == "yes") shuffle = true;
				else if (opt_arg == "no") shuffle = false;
				else if (proc_id == 0) std::cerr << "warning: shuffle option expects 'yes' or 'no'\n";
				break;

			// set vectorized input to resume hierarchy run
			case 'R':
				resume_input_fname = opt_arg;
				resume = true;
				break;

			// set output directory
			case 'o':
				output_dir = opt_arg;
				break;

			// set output format
			case 'f':
				if (opt_arg == "clusters") output_format = hof_clusters;
				else if (opt_arg == "vectorized") output_format = hof_vectorized;
				else if (opt_arg == "none") output_format = hof_none;
				else if (proc_id == 0) std::cerr << "warning: unrecognized format descriptor: " << opt_arg << "\n";
				break;

			// print usage and exit
			case 'h':
				if (proc_id == 0)
				{
#ifdef WITH_MPI
					//std::cout << "usage: mpirun -n <numprocs> HMACTool <inputfile> -s <sigmas> [-r <minpart>]\n\n";
					opt_handler.printUsage("usage: mpirun [-n numprocs] HMACTool -s sigmas [-r minpart] inputfile");
#else
					//std::cout << "usage: HMACTool <inputfile> -s <sigmas> [-r <minpart>]\n\n";
					opt_handler.printUsage("usage: HMACTool -s sigmas [-r minpart] inputfile");
#endif
/*
					std::cout << "command line options:\n";
					std::cout << "\t--help, -h : print this usage dialog and exit\n";
					std::cout << "\t--verbose, -v : print verbose runtime information to stdout\n";
					std::cout << "\t--quiet, -q : suppress output to stdout\n";
					std::cout << " *REQ*\t--sigmas, -s <sigmafile> : specify path to file containing sigmas (no default, required argument)\n";
					std::cout << "\t--recursive, -r <minpart> : perform recursive HMAC using minimum partition size specified by argument (default no recursion)\n";
					std::cout << "\t--normalize, -N [(yes)/no] : perform normalization of data prior to computation (default yes)\n";
					std::cout << "\t--shuffle, -S [(yes)/no] : perform shuffling of data prior to partitioning (default yes)\n";
					std::cout << "\t--resume, -R <inputfile> : specify vectorized output to resume hierarchy construction\n";
					std::cout << "\t--output-dir, -o <outputdir> : set directory to write outputs for each sigma (default ./)\n";
					std::cout << "\t--output-format, -f <outputfmt> : set format of output files [(clusters), vectorized, none]\n";
					std::cout << "\n";

					std::cout << "\t<required argument>   [optional argument]   (default argument)\n\n";
*/
				}

			// exit with unrecognized option
			case '?':
				hmac_exit(0);
				break;
		};
	};

	// initialize timer
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	// no input file specified
	//if (optind == argc)
	if (!opt_handler.getNonOptArg(input_fname))
	{
		if (proc_id == 0) std::cerr << "error: no input file specified\n";
		hmac_exit(1);
	}

	// multiple input files specified
	//else if (optind + 1 < argc)
	std::string dummy_name;
	if (opt_handler.getNonOptArg(dummy_name))
	{
		if (proc_id == 0)
		{
			std::cerr << "error: multiple input files specified: '" << input_fname << "'";
			do { std::cerr << " '" << dummy_name << "'"; }
			while (opt_handler.getNonOptArg(dummy_name));
			std::cerr << "\n";
		}
/*
		if (proc_id == 0)
		{
			std::cerr << "error: multiple input files specified:";
			while (optind < argc) std::cerr << " '" << argv[optind++] << "'";
			std::cerr << "\n";
		}
*/
		hmac_exit(2);
	}

	//// get input file name
	//else input_fname = argv[optind];

	// no sigma specified
	if (sigmas.empty())
	{
		if (proc_id == 0) std::cerr << "error: no sigmas specified\n";
		hmac_exit(3);
	}

#ifdef WITH_BOOST
	/*// output directory does not exist
	if (!boost::filesystem::is_directory(output_dir))
	{
		if (proc_id == 0) std::cerr << "error: output directory does not exist: " << output_dir << "\n";
		hmac_exit(7);
	}*/

	// NOTE: (attempted) quick fix to checking permissions and creating output directory as necessary

	// output directory does not exist
	if (!boost::filesystem::is_directory(output_dir))
	{
#ifdef WITH_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		// attempt to create directory
		if (proc_id == 0) boost::filesystem::create_directories(output_dir);

#ifdef WITH_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		if (!boost::filesystem::is_directory(output_dir))
		{
			if (proc_id == 0) std::cerr << "error: could not create output directory: " << output_dir << "\n";
			hmac_exit(7);
		}
	}

	// check permissions of existing directory
	else
	{
		// attempt to write to dummy file in directory
		std::string dummy_fname = output_dir + "/." + std::to_string((long long) proc_id);
		std::ofstream dummy_file(dummy_fname);
		dummy_file << "hello\n";

		// existing directory is not writable
		if (dummy_file.fail())
		{
			if (proc_id == 0) std::cerr << "error: output directory is not writable: " << output_dir << "\n";
			hmac_exit(9);
		}

		dummy_file.close();
		boost::filesystem::remove(dummy_fname);
	}
#endif

	// open dataset
	std::vector<double> coords;
	unsigned long num_coords, num_fields;
	if (!open_dataset(input_fname, coords, num_coords, num_fields))
		hmac_exit(4);

	// print runtime information
	if (verbose)
	{
		if (proc_id == 0)
		{
			std::cout << "input file name = " << input_fname << "\n";
			std::cout << "num coords = " << num_coords << "\n";
			std::cout << "num fields = " << num_fields << "\n";
			std::list<double>::iterator sigma = sigmas.begin();
			if (sigmas.size() == 1) std::cout << "sigma = " << *sigma;
			else std::cout << "sigmas = " << *sigma;
			while (++sigma != sigmas.end()) std::cout << ", " << *sigma;
			std::cout << "\n";
			std::cout << "recursion = ";
			if (min_part_size == 0) std::cout << "false\n";
			else std::cout << "true (" << min_part_size << ")\n";
			std::cout << "data normalization = " << (normalize ? "true" : "false") << "\n";
			std::cout << "data shuffling = " << (shuffle ? "true" : "false") << "\n";
			std::cout << "output directory = " << ((output_dir == ".") ? "(local)" : output_dir) << "\n";
		}
	}

	// initialize index vector
	std::vector<unsigned long> shuffle_indices(num_coords);
	for (unsigned long i = 0; i < num_coords; ++i)
		shuffle_indices[i] = i;

	// perform data shuffling
	std::list<unsigned long> shuffle_moves;
	if (shuffle)
	{
		for (unsigned long i = 0; i < num_coords; ++i)
		{
			unsigned long j = rand() % num_coords;
			shuffle_moves.push_back(j);
			if (i == j) continue;

			// swap indices
			std::swap(shuffle_indices[i], shuffle_indices[j]);

			std::vector<double> tmp(num_fields);
			std::vector<double>::iterator coord_i = coords.begin() + num_fields * i;
			std::vector<double>::iterator coord_j = coords.begin() + num_fields * j;
			std::vector<double>::iterator coord_tmp = tmp.begin();
 
			// swap coordinates
			std::copy_n(coord_i, num_fields, coord_tmp);
			std::copy_n(coord_j, num_fields, coord_i);
			std::copy_n(coord_tmp, num_fields, coord_j);
		}
	}

	// perform data normalization
	std::vector<double> scalars(num_fields, 1.0), averages(num_fields, 0.0);
	if (normalize)
	{
		// calculate averages for each field
		for (unsigned long i = 0; i < num_coords; ++i)
			for (unsigned long j = 0; j < num_fields; ++j)
				averages[j] += coords[num_fields * i + j];
		for (unsigned long i = 0; i < num_fields; ++i)
			averages[i] /= num_coords;

		if (verbose && proc_id == 0)
		{
			std::cout << "averages =";
			for (unsigned long i = 0; i < num_fields; ++i)
				std::cout << " " << averages[i];
			std::cout << "\n";
		}

		// calculate variance
		std::vector<double> variance(num_fields);
		for (unsigned long i = 0; i < num_fields; ++i)
			variance[i] = 0.0;
		for (unsigned long i = 0; i < num_coords; ++i)
			for (unsigned long j = 0; j < num_fields; ++j)
				variance[j] += (coords[num_fields * i + j] - averages[j]) * (coords[num_fields * i + j] - averages[j]);

		// calculate standard deviation
		for (unsigned long i = 0; i < num_fields; ++i)
			scalars[i] = sqrt(variance[i] / num_coords);

		if (verbose && proc_id == 0)
		{
			std::cout << "scalars =";
			for (unsigned long i = 0; i < num_fields; ++i)
				std::cout << " " << scalars[i];
			std::cout << "\n";
		}

		// normalize data
		for (unsigned long i = 0; i < num_coords; ++i)
			for (unsigned long j = 0; j < num_fields; ++j)
				coords[num_fields * i + j] = (coords[num_fields * i + j] - averages[j]) / scalars[j];
	}

	// initialize hmac dataset struct
	hmac_dataset dataset;
	dataset.begin = coords.begin();
	dataset.end = coords.end();
	dataset.num_coords = num_coords;
	dataset.num_fields = num_fields;

	// NOTE: resume option skips the initial hierarchy level
	if (!resume) {

		// initialize process workload
		unsigned long num_to_cluster = num_coords / num_procs;
		unsigned long start_index = num_to_cluster * proc_id;
		if (proc_id < num_coords % num_procs) ++num_to_cluster;
		start_index += std::min(proc_id, (int) (num_coords % num_procs));

		// initialize to_cluster dataset
		hmac_dataset to_cluster;
		to_cluster.begin = dataset.begin + num_fields * start_index;
		to_cluster.end = to_cluster.begin + num_fields * num_to_cluster;
		to_cluster.num_coords = num_to_cluster;
		to_cluster.num_fields = num_fields;

		// print preprocessing time
		std::chrono::system_clock::time_point preprocess;
		if (verbose && proc_id == 0)
		{
			preprocess = std::chrono::system_clock::now();
			std::chrono::duration<double> time_diff = preprocess - start;
			std::cout << "pre-processing time = " << time_diff.count() << " seconds\n";
		}

		std::map<hmac_mode, std::list<unsigned long>, mode_lt> clusters;

		// run non-recursive MAC on to_cluster with full dataset
		if (min_part_size >= num_coords || min_part_size < 1)
		{
			if (verbose && proc_id == 0) std::cout << "running non-recursive MAC with full dataset\n\n";
			ModalAC(dataset, to_cluster, sigmas.front(), clusters);
		}

		else
		{
			std::map<hmac_mode, std::list<unsigned long>, mode_lt> tmp_clusters;

			// run non-recursive MAC on to_cluster with windowed dataset
			if (min_part_size >= num_to_cluster)
			{
				// create data window centered around to_cluster
				signed long left = start_index - (min_part_size - num_to_cluster) / 2;
				signed long right = left + min_part_size;

				// account for left margin
				if (left < 0)
				{
					right -= left;
					left = 0;
				}

				// account for right margin
				else if (right > num_coords)
				{
					left -= (right - num_coords);
					right = num_coords;
				}

				// initialize window
				hmac_dataset window;
				window.begin = dataset.begin + num_fields * left;
				window.end = dataset.begin + num_fields * right;
				window.num_coords = min_part_size;
				window.num_fields = num_fields;

				// run MAC on to_cluster with windowed dataset
				if (verbose && proc_id == 0) std::cout << "running non-recursive MAC with windowed dataset\n\n";
				ModalAC(window, to_cluster, sigmas.front(), tmp_clusters);
			}

			// run recursive MAC on to_cluster partition
			else
			{
				if (verbose && proc_id == 0) std::cout << "running recursive MAC\n\n";
				RecursiveMAC(to_cluster, sigmas.front(), min_part_size, tmp_clusters);
			}

			// run MEM on modes using full dataset
			std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator tmp_cluster;
			for (tmp_cluster = tmp_clusters.begin(); tmp_cluster != tmp_clusters.end(); ++tmp_cluster)
			{
				std::vector<double> start(tmp_cluster->first), mode(dataset.num_fields);
				ModalEM(dataset, sigmas.front(), start, mode);
				ClusterMapMerge(clusters, mode, tmp_cluster->second);
			}
		}

		// print level 1 processing time
		std::chrono::system_clock::time_point postprocess;
		if (verbose && proc_id == 0)
		{
			postprocess = std::chrono::system_clock::now();
			std::chrono::duration<double> time_diff = postprocess - preprocess;
			std::cout << "level 1 processing time = " << time_diff.count() << " seconds\n";
		}

		// NOTE: every proc at this point has fully clustered partition in "clusters"

#ifdef WITH_MPI
		unsigned long index_offset = num_to_cluster;

		// receive results from "child" process
		for (int i = 1; i < num_procs; i <<= 1)
		{
			int child_id = proc_id | i;
			if (child_id == proc_id || child_id >= num_procs) break;

			// receive vectorized results
			std::vector<double> tmp_modes;
			std::vector<unsigned long> tmp_indices;
			RecvModesAndIndices(child_id, tmp_modes, tmp_indices);

			// transform to cluster map
			std::map<hmac_mode, std::list<unsigned long>, mode_lt> tmp_clusters;
			ModesAndIndices2Clusters(tmp_modes, tmp_indices, num_fields, tmp_clusters);

			// insert new clusters or merge with existing clusters
			std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator tmp_cluster;
			for (tmp_cluster = tmp_clusters.begin(); tmp_cluster != tmp_clusters.end(); ++tmp_cluster)
			{
				// offset index ranges of cluster
				std::list<unsigned long>::iterator index;
				for (index = tmp_cluster->second.begin(); index != tmp_cluster->second.end(); ++index)
					*index += index_offset;

				ClusterMapMerge(clusters, tmp_cluster->first, tmp_cluster->second);
			}

			index_offset += tmp_indices.size();
		}

		// send results to "parent" process
		if (proc_id > 0)
		{
			// get least significant bit of proc_id
			int lsb = 1;
			while (!(lsb & proc_id)) lsb <<= 1;

			int parent_id = proc_id - lsb;

			// vectorize results for MPI transfer
			std::vector<double> tmp_modes;
			std::vector<unsigned long> tmp_indices;
			Clusters2ModesAndIndices(clusters, index_offset, tmp_modes, tmp_indices);

			// send results to parent
			SendModesAndIndices(parent_id, tmp_modes, tmp_indices);
		}
#endif

		// NOTE: top process now contains fully clustered dataset in "clusters"

		// print merge time
		std::chrono::system_clock::time_point merge;
		if (verbose && proc_id == 0)
		{
			merge = std::chrono::system_clock::now();
			std::chrono::duration<double> time_diff = merge - postprocess;
			std::cout << "level 1 merge time = " << time_diff.count() << " seconds\n";
		}

		/*// print details of level to stdout
		if (!quiet && proc_id == 0)
		{
			std::cout << std::endl << clusters.size() << " modes on level 1 (sigma = " << sigmas.front() << ")";

			// don't print modes for levels of 10 or more clusters
			if (clusters.size() > 9) std::cout << "\n";

			// print modes for reasonably small number of clusters
			else
			{
				std::cout << ":\n";

				std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator cluster;
				for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster)
				{
					// print mode
					for (unsigned long i = 0; i < num_fields; ++i)
						std::cout << scalars[i] * cluster->first[i] + averages[i] << " ";

					// print membership size
					std::cout << "(" << cluster->second.size() << ")\n";
				}
			}

			std::cout << std::endl;
		}*/

		resume_input_fname = output_dir + "/level-1.txt";
		if (verbose && proc_id == 0) std::cout << "writing level 1 result to '" << resume_input_fname << "'...\n";

		// write output files
		if (proc_id == 0)
		{
			output_format = hof_vectorized; // TODO create clusters-based output
			if (output_format == hof_vectorized)
			{
				std::vector<double> modes;
				std::vector<unsigned long> indices;
				Clusters2ModesAndIndices(clusters, num_coords, modes, indices);

				std::ofstream fout(output_dir + "/level-1.txt");

				// print header
				fout << clusters.size() << " " << num_coords << " " << num_fields << "\n";

				// print un-normalized modes
				for (unsigned long i = 0; i < clusters.size(); ++i)
				{
					for (unsigned long j = 0; j < num_fields; ++j)
						fout << scalars[j] * modes[num_fields * i + j] + averages[j] << " ";

					fout << "\n";
				}

				// map shuffled indices to correct locations
				std::vector<unsigned long> real_indices(num_coords);
				for (unsigned long i = 0; i < num_coords; ++i)
					real_indices[shuffle_indices[i]] = indices[i];

				// print indices
				for (unsigned long i = 0; i < num_coords; ++i)
					fout << real_indices[i] << "\n";
			}
		}

		// print details of level to stdout
		if (!quiet && proc_id == 0)
		{
			std::cout << std::endl << clusters.size() << " modes on level 1 (sigma = " << sigmas.front() << ")";

			// don't print modes for levels of 10 or more clusters
			if (clusters.size() > 9) std::cout << "\n";

			// print modes for reasonably small number of clusters
			else
			{
				std::cout << ":\n";

				std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator cluster;
				for (cluster = clusters.begin(); cluster != clusters.end(); ++cluster)
				{
					// print mode
					for (unsigned long i = 0; i < num_fields; ++i)
						std::cout << scalars[i] * cluster->first[i] + averages[i] << " ";

					// print membership size
					std::cout << "(" << cluster->second.size() << ")\n";
				}
			}

			std::cout << std::endl;
		}

		sigmas.pop_front();

	// NOTE: end of resume clause
	}

	int level;

	// begin hierarchy construction with level 1 result
	if (!resume)
	{
		level = 1;
		//resume_input_fname = output_dir + "/level-1.txt";
	}

	// extract level number from file name
	else
	{
		// get file name without directory
		std::string tmp;
		size_t fname_pos = resume_input_fname.find_last_of("/") + 1;
		if (fname_pos == resume_input_fname.size())
			tmp = resume_input_fname;
		else tmp = resume_input_fname.substr(fname_pos);

		// check input name
		if (tmp.size() < 11)
		{
			if (proc_id == 0) std::cerr << "error: resume file does not follow naming convention\n";
			hmac_exit(5);
		}

		// level-(value).txt
		level = (int) std::stoul(tmp.substr(6, tmp.size() - 10));
	}

	while (!sigmas.empty())
	{
#ifdef WITH_MPI
		// make sure process 0 completes previous file writes
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		//resume_input_fname = "level-" + std::to_string(level++) + ".txt";
		++level;

		// open vectorized previous result
		unsigned long num_modes;
		std::vector<double> modes;
		std::vector<unsigned long> indices;
		if (!open_vectorized(resume_input_fname, modes, indices, num_modes, num_coords, num_fields))
			hmac_exit(6);

		// perform data normalization
		if (normalize)
			for (unsigned long i = 0; i < num_modes; ++i)
				for (unsigned long j = 0; j < num_fields; ++j)
					modes[num_fields * i + j] = (modes[num_fields * i + j] - averages[j]) / scalars[j];

		// reinitialize to_cluster dataset
		unsigned long num_to_cluster = num_modes / num_procs;
		unsigned long start_index = num_to_cluster * proc_id;
		if (proc_id < num_modes % num_procs) ++num_to_cluster;
		start_index += std::min(proc_id, (int) (num_modes % num_procs));

		// reinitialize to_cluster dataset
		hmac_dataset to_cluster;
		to_cluster.begin = modes.begin() + num_fields * start_index;
		to_cluster.end = to_cluster.begin + num_fields * num_to_cluster;
		to_cluster.num_coords = num_to_cluster;
		to_cluster.num_fields = num_fields;

		std::chrono::system_clock::time_point preprocess = std::chrono::system_clock::now();

		// perform Modal AC on previous modes
		std::map<hmac_mode, std::list<unsigned long>, mode_lt> clusters;
		ModalAC(dataset, to_cluster, sigmas.front(), clusters);

		// print level processing time
		std::chrono::system_clock::time_point postprocess;
		if (verbose && proc_id == 0)
		{
			postprocess = std::chrono::system_clock::now();
			std::chrono::duration<double> time_diff = postprocess - preprocess;
			std::cout << "level " << level << " processing time = " << time_diff.count() << " seconds\n";
		}

#ifdef WITH_MPI
		unsigned long index_offset = num_to_cluster;

		// receive results from "child" process
		for (int i = 1; i < num_procs; i <<= 1)
		{
			int child_id = proc_id | i;
			if (child_id == proc_id || child_id >= num_procs) break;

			// receive vectorized results
			std::vector<double> tmp_modes;
			std::vector<unsigned long> tmp_indices;
			RecvModesAndIndices(child_id, tmp_modes, tmp_indices);

			// transform to cluster map
			std::map<hmac_mode, std::list<unsigned long>, mode_lt> tmp_clusters;
			ModesAndIndices2Clusters(tmp_modes, tmp_indices, num_fields, tmp_clusters);

			// insert new clusters or merge with existing clusters
			std::map<hmac_mode, std::list<unsigned long>, mode_lt>::iterator tmp_cluster;
			for (tmp_cluster = tmp_clusters.begin(); tmp_cluster != tmp_clusters.end(); ++tmp_cluster)
			{
				// offset index ranges of cluster
				std::list<unsigned long>::iterator index;
				for (index = tmp_cluster->second.begin(); index != tmp_cluster->second.end(); ++index)
					*index += index_offset;

				ClusterMapMerge(clusters, tmp_cluster->first, tmp_cluster->second);
			}

			index_offset += tmp_indices.size();
		}

		// send results to "parent" process
		if (proc_id > 0)
		{
			// get least significant bit of proc_id
			int lsb = 1;
			while (!(lsb & proc_id)) lsb <<= 1;

			int parent_id = proc_id - lsb;

			// vectorize results for MPI transfer
			std::vector<double> tmp_modes;
			std::vector<unsigned long> tmp_indices;
			Clusters2ModesAndIndices(clusters, index_offset, tmp_modes, tmp_indices);

			// send results to parent
			SendModesAndIndices(parent_id, tmp_modes, tmp_indices);
		}

		// print merge time
		std::chrono::system_clock::time_point merge;
		if (verbose && proc_id == 0)
		{
			merge = std::chrono::system_clock::now();
			std::chrono::duration<double> time_diff = merge - postprocess;
			std::cout << "level " << level << " merge time = " << time_diff.count() << " seconds\n";
		}
#endif

		// update modes and indices
		std::vector<unsigned long> membership_size(clusters.size(), 0);
		if (proc_id == 0)
		{
			std::vector<double> mode_modes;
			std::vector<unsigned long> mode_indices;
			Clusters2ModesAndIndices(clusters, num_modes, mode_modes, mode_indices);

			// map old indices to new indices
			for (unsigned long i = 0; i < num_coords; ++i)
			{
				indices[i] = mode_indices[indices[i]];
				membership_size[indices[i]]++;
			}

			// update modes vector
			modes.swap(mode_modes);
		}

		/*// print details of level to stdout
		if (!quiet && proc_id == 0)
		{
			std::cout << std::endl << clusters.size() << " modes on level " << level << " (sigma = " << sigmas.front() << ")";

			// don't print modes for levels of 10 or more clusters
			if (clusters.size() > 9) std::cout << "\n";

			// print modes for reasonably small number of clusters
			else
			{
				std::cout << ":\n";

				for (unsigned long i = 0; i < clusters.size(); ++i)
				{
					// print mode
					for (unsigned long j = 0; j < num_fields; ++j)
						std::cout << scalars[j] * modes[num_fields * i + j] + averages[j] << " ";

					// print membership size
					std::cout << "(" << membership_size[i] << ")\n";
				}
			}

			std::cout << std::endl;
		}*/

		resume_input_fname = output_dir + "/level-" + std::to_string((long long) level) + ".txt";
		if (verbose && proc_id == 0) std::cout << "writing level " << level << " result to '" << resume_input_fname << "'...\n";

		// write output files
		if (proc_id == 0)
		{
			output_format = hof_vectorized; // TODO remove
			if (output_format == hof_vectorized)
			{
				std::ofstream fout(resume_input_fname);

				// print header
				fout << clusters.size() << " " << num_coords << " " << num_fields << "\n";

				// print un-normalized modes
				for (unsigned long i = 0; i < clusters.size(); ++i)
				{
					for (unsigned long j = 0; j < num_fields; ++j)
						fout << scalars[j] * modes[num_fields * i + j] + averages[j] << " ";

					fout << "\n";
				}

				// print indices
				for (unsigned long i = 0; i < num_coords; ++i)
					fout << indices[i] << "\n";
			}
		}

		// print details of level to stdout
		if (!quiet && proc_id == 0)
		{
			std::cout << std::endl << clusters.size() << " modes on level " << level << " (sigma = " << sigmas.front() << ")";

			// don't print modes for levels of 10 or more clusters
			if (clusters.size() > 9) std::cout << "\n";

			// print modes for reasonably small number of clusters
			else
			{
				std::cout << ":\n";

				for (unsigned long i = 0; i < clusters.size(); ++i)
				{
					// print mode
					for (unsigned long j = 0; j < num_fields; ++j)
						std::cout << scalars[j] * modes[num_fields * i + j] + averages[j] << " ";

					// print membership size
					std::cout << "(" << membership_size[i] << ")\n";
				}
			}

			std::cout << std::endl;
		}

		sigmas.pop_front();
	}

	// print total run time
	if ((verbose || !quiet) && proc_id == 0)
	{
		std::chrono::system_clock::time_point finish = std::chrono::system_clock::now();
		std::chrono::duration<double> time_diff = finish - start;
		std::cout << "total run time = " << time_diff.count() << " seconds\n";
	}

	hmac_exit(0);
}

bool open_dataset(const std::string &file_name, std::vector<double> &coords, unsigned long &num_coords, unsigned long &num_fields)
{
	coords.clear();
	num_coords = num_fields = 0;
	std::ifstream fin(file_name);

	// unable to open file
	if (!fin.is_open())
	{
		if (proc_id == 0) std::cerr << "error: could not open file: " << file_name << std::endl;
		return false;
	}

	// get one coord per line
	std::string line;
	while (std::getline(fin, line))
	{
		unsigned long nf = 0;
		std::stringstream ss(line);

		while (!ss.eof())
		{
			char c = ss.peek();
			if (isdigit(c) || c == '-' || c == '.')
			{
				double val;
				ss >> val;
				coords.push_back(val);
				++nf;
			}

			else ss.get();
		}

		// empty line
		if (nf == 0) continue;

		++num_coords;
		
		// initialize num fields
		if (num_fields == 0) num_fields = nf;

		// conflicting number of fields
		else if (num_fields != nf)
		{
			if (proc_id == 0) std::cerr << "error: conflicting number of fields on coord " << num_coords << std::endl;
			return false;
		}
	}

	return true;
}

bool open_vectorized(const std::string &file_name, std::vector<double> &modes, std::vector<unsigned long> &indices, unsigned long &num_modes, unsigned long &num_coords, unsigned long &num_fields)
{
	modes.clear();
	indices.clear();
	std::ifstream fin(file_name);

	// unable to open file
	if (!fin.is_open())
	{
		if (proc_id == 0) std::cerr << "error: could not open file: " << file_name << std::endl;
		return false;
	}

	// get header values
	fin >> num_modes;
	fin >> num_coords;
	fin >> num_fields;

	modes.resize(num_fields * num_modes);
	indices.resize(num_coords);

	// get modes
	for (unsigned long i = 0; i < num_modes; ++i)
		for (unsigned long j = 0; j < num_fields; ++j)
			fin >> modes[num_fields * i + j];

	// get indices
	for (unsigned long i = 0; i < num_coords; ++i)
		fin >> indices[i];

	return true;
}

#ifdef WITH_MPI
void SendModesAndIndices(int to_proc_id, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	// send buffer sizes
	int transfer_sizes[] = { (int) modes.size(), (int) indices.size() };
	MPI_Send(transfer_sizes, 2, MPI_INT, to_proc_id, 1, MPI_COMM_WORLD);

	// send buffers
	MPI_Send(modes.data(), transfer_sizes[0], MPI_DOUBLE, to_proc_id, 2, MPI_COMM_WORLD);
	MPI_Send(indices.data(), transfer_sizes[1], MPI_UNSIGNED_LONG, to_proc_id, 3, MPI_COMM_WORLD);
}

void RecvModesAndIndices(int from_proc_id, std::vector<double> &modes, std::vector<unsigned long> &indices)
{
	// receive buffer sizes
	int transfer_sizes[2];
	MPI_Recv(transfer_sizes, 2, MPI_INT, from_proc_id, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// receive buffers
	modes.resize(transfer_sizes[0]);
	indices.resize(transfer_sizes[1]);
	MPI_Recv(modes.data(), transfer_sizes[0], MPI_DOUBLE, from_proc_id, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(indices.data(), transfer_sizes[1], MPI_UNSIGNED_LONG, from_proc_id, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
#endif

