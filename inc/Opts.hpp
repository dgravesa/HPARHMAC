/*
	Daniel Graves
	Opts.hpp
*/

#ifndef OPTS_HPP
#define OPTS_HPP

#include <string>
#include <iostream>
#include <list>
#include <unordered_map>
#include <cctype>
#include <cstring>
#include <algorithm>

enum OptType { REQUIRED_ARG, OPTIONAL_ARG, NO_ARG };

struct Option
{
	const char *long_opt;
	char short_opt;
	int type;
	const char *description;
};

class OptHandler
{
public:
	OptHandler(const Option *options, int argc, char **argv) { init(options, argc, argv); }
	~OptHandler() { ; }

	// initialize option handler, same arguments as constructor
	void init(const Option *options, int argc, char **argv);

	// get next option and argument, returns false if no more options are available, true otherwise
	bool getOpt(char &opt_code, std::string &opt_arg);

	// get next non-option argument, returns false if no more non-option arguments are available, true otherwise
	bool getNonOptArg(std::string &non_opt_arg);

	// print usage dialog based on provided options and arguments
	void printUsage(const std::string &usage_line = "") const;

private:
	const Option *opts;
	std::unordered_map<char, const Option *> short_opts;
	std::unordered_map<std::string, const Option *> long_opts;
	std::list<std::string> non_opt_args;

	int arg_index, arg_count;
	char **args, *arg_ptr;
	bool new_arg;
};

void OptHandler::init(const Option *options, int argc, char **argv)
{
	opts = options;
	arg_count = argc;
	args = argv;
	arg_index = 1;
	new_arg = true;

	// no options provided
	if (opts == NULL)
	{
		// add all arguments to non-option arguments list
		while (arg_index < arg_count)
			non_opt_args.push_back(args[arg_index++]);
	}

	else
	{
		// initialize option hash tables
		for (int i = 0; opts[i].short_opt != 0; ++i)
		{
			if (std::isalnum(opts[i].short_opt))
				short_opts.insert(std::pair<char, const Option *>(opts[i].short_opt, opts + i));

			if (opts[i].long_opt != NULL)
				long_opts.insert(std::pair<std::string, const Option *>(opts[i].long_opt, opts + i));
		}
	}
}

bool OptHandler::getOpt(char &opt_code, std::string &opt_arg)
{
	opt_code = '?';
	opt_arg = "";

	if (arg_index >= arg_count)
		return false;

	if (new_arg)
	{
		arg_ptr = args[arg_index];

		// argument is an option
		if (arg_ptr[0] == '-' && arg_ptr[1] != 0)
		{
			if (arg_ptr[1] == '-')
			{
				// stop parsing options
				if (arg_ptr[2] == 0)
				{
					while (++arg_index < arg_count)
						non_opt_args.push_back(args[arg_index]);

					return false;
				}

				// new long option argument
				else
				{
					arg_ptr += 2;

					// get possible optional argument
					char *split = strchr(arg_ptr, '=');
					if (split != NULL)
					{
						*split = 0;
						opt_arg = split + 1;
					}

					std::unordered_map<std::string, const Option *>::iterator opt_find = long_opts.find(arg_ptr);

					// long option not recognized
					if (opt_find == long_opts.end())
					{
						std::cerr << "error: unrecognized long option '" << arg_ptr << "'\n";
						arg_index = arg_count;
						return true;
					}

					// unexpected option argument error
					if (opt_find->second->type == NO_ARG && !opt_arg.empty())
					{
						std::cerr << "error: long option '" << arg_ptr << "' expects no argument\n";
						arg_index = arg_count;
						return true;
					}

					if (opt_find->second->type == REQUIRED_ARG && opt_arg.empty())
					{
						// expected argument error
						if (++arg_index >= arg_count)
						{
							std::cerr << "error: long option '" << arg_ptr << "' requires an argument\n";
							return true;
						}

						// set option argument as next argument
						opt_arg = args[arg_index];
					}

					opt_code = opt_find->second->short_opt;

					++arg_index;
					return true;
				}
			}

			// new short options argument
			else
			{
				++arg_ptr;
				new_arg = false;
				return getOpt(opt_code, opt_arg);
			}
		}

		// argument is not an option
		else
		{
			non_opt_args.push_back(args[arg_index++]);
			return getOpt(opt_code, opt_arg);
		}
	}

	// get additional short options from argument
	else
	{
		std::unordered_map<char, const Option *>::iterator opt_find = short_opts.find(*arg_ptr);

		// short option not recognized
		if (opt_find == short_opts.end())
		{
			std::cerr << "error: short option '" << *arg_ptr << "' not recognized\n";
			arg_index = arg_count;
			return true;
		}

		// set option argument as remainder of argument
		if (opt_find->second->type != NO_ARG)
			opt_arg = arg_ptr + 1;

		// set option argument as next argument
		if (opt_find->second->type == REQUIRED_ARG && opt_arg.empty())
		{
			// expected argument error
			if (++arg_index >= arg_count)
			{
				std::cerr << "error: short option '" << *arg_ptr << "' expects an argument\n";
				return true;
			}

			opt_arg = args[arg_index];
		}

		++arg_ptr;

		// move to next argument
		if (opt_find->second->type != NO_ARG || *arg_ptr == 0)
		{
			++arg_index;
			new_arg = true;
		}

		opt_code = opt_find->second->short_opt;

		return true;
	}
}

bool OptHandler::getNonOptArg(std::string &non_opt_arg)
{
	non_opt_arg = "";

	if (non_opt_args.empty())
		return false;

	non_opt_arg = non_opt_args.front();
	non_opt_args.pop_front();

	return true;
}

void OptHandler::printUsage(const std::string &usage_str) const
{
	// print usage line
	if (!usage_str.empty())
		std::cout << usage_str << "\n";

	if (opts == NULL) return;

	// get maximum long option name for fill
	size_t fill_width = 0;
	for (int i = 0; opts[i].short_opt != 0; ++i)
	{
		size_t this_fill = 4;

		if (std::isalnum(opts[i].short_opt))
		{
			this_fill += 2;
			if (opts[i].type != NO_ARG)
				this_fill += 7;
			if (opts[i].long_opt != NULL)
				this_fill += 2;
		}

		if (opts[i].long_opt != NULL)
		{
			this_fill += 2 + strlen(opts[i].long_opt);

			if (opts[i].type != NO_ARG)
				this_fill += 8;
		}

		fill_width = std::max(fill_width, this_fill);
	}
	//fill_width += 10;

	// print command line options
	std::cout << "\ncommand line options:\n";
	for (int i = 0; opts[i].short_opt != 0; ++i)
	{
		std::string opt_name = "  ";

		// add short option to name string
		if (std::isalnum(opts[i].short_opt))
		{
			opt_name += "-";
			opt_name += opts[i].short_opt;

			switch (opts[i].type)
			{
				case REQUIRED_ARG:
					opt_name += " VALUE";
					break;

				case OPTIONAL_ARG:
					opt_name += "[VALUE]";
					break;

				default:
					break;
			};

			if (opts[i].long_opt != NULL)
				opt_name += ", ";
		}

		// add long option to name string
		if (opts[i].long_opt != NULL)
		{
			opt_name += "--";
			opt_name += opts[i].long_opt;

			switch (opts[i].type)
			{
				case REQUIRED_ARG:
					opt_name += " VALUE";
					break;

				case OPTIONAL_ARG:
					opt_name += "[=VALUE]";
					break;

				default:
					break;
			};
		}

		// fill width
		size_t fill = fill_width - opt_name.size();
		opt_name.append(fill, ' ');
		std::cout << opt_name;

		// print description of option
		if (opts[i].description != NULL)
			std::cout << ": " << opts[i].description;

		std::cout << "\n";
	}

	std::cout << "\n";
}

#endif

