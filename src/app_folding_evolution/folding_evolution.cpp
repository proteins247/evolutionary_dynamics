/* 
 * folding_evolution version b0.0.14-draft4 (drafting)
 *
 * branch: stability_only
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <csignal> 		// eventually need to think about handling interrupts
#include <cstring>
#include <cstdio>
#include <cmath>
#include <map>
#include <limits>

#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <dirent.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>


#include <hdf5.h>
#include <mpi.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#define CPLUSPLUS __cplusplus
#undef __cplusplus
// When compiling Random123 header files for rng.h, we don't need the
//   C++ features (and they're not compatible with extern "C"), so we
//   temporarily undefine __cplusplus, saving its value in CPLUSPLUS.
//   Alternatively, we could build rng.o and gencode.o using g++, and
//   then we wouldn't have to wrap with extern "C".

extern "C"
{
#include "../gencode.h"
#include "../rng.h"
}

#define __cplusplus CPLUSPLUS
#undef CPLUSPLUS


static const std::string helptext =
    "HELP:\n"
    "folding_evolution\n"
    "Requires latPack and MPI\n"
    "\n"
    "Usage: mpirun -np N folding_evolution [OPTIONS] SEQUENCE N_GENS \\\n"
    "         --native-fold NFOLD --speed-params FILE\n"
    "\n"
    "Carry out folding evolutionary simulations. This is a parallel program\n"
    "that carries out as many folding simulations each generation as there\n"
    "are available processors.\n"
    "\n"
    "Folding simulations are used to evaluate the fitness of a a genetic\n"
    "sequence. Because such a fitness measurement is noisy, sequences may\n"
    "show fitnesses higher or lower than the actual values. Hence with each\n"
    "generation, the current sequence is also reevaluated.\n"
    "\n"
    "  -n, --native-fold=CONF    REQUIRED. The native conformation of the\n"
    "                            as a latPack-style move sequence.\n"
    "  -a, --activity-constant   Set constant in fitness function.\n"
    "                            Default=0.25\n"
    "  -d, --debug               Allow error messages from all processors.\n"
    "  -f, --from-checkpoint     Resume simulation from checkpoint file.\n"
    "  -k, --degradation-param=K Value setting timescale for protein\n"
    "                            degradation (penalizes slow folding and\n"
    "                            low stability).\n"
    "  -m, --mutation-mode=M     One of 0, 1, or 2, indicating restriction to\n"
    "                            synonymous mutation, nonsynonymous mutation,\n"
    "                            or allow both kinds of mutations (default).\n"
    "      --new-ngens=N         When resuming from checkpoint, change the\n"
    "                            target number of generations."    
    "  -o, --out-path=FILE       Output will be put in this directory.\n"
    "                            Default=./out\n"
    "  -p, --population-size=N   Size of population for evolutionary\n"
    "                            dynamics. Default=500\n"
    "      --random-codons       If SEQUENCE is of amino acids, corresponding\n"
    "                            RNA sequence codons will be randomly chosen.\n"
    "                            (Default behavior is to use fastest codons.)\n"
    "  -r, --seed=N              RNG seed. Default=1\n"
    "  -t, --temperature=T       Temperature of latfold simulations.\n"
    "                            Default=0.2\n"
    "      --help   Display this help and exit.\n"
    "\n"
    "Format of the output files:\n"
    "\nThe primary output is a simulation log that records the generation\n"
    "number, protein sequence, nucleic acid sequence, fitness, and whether\n"
    "the mutation was accepted. Trajectory files are also saved.\n"
    "\n"
    ;

enum MutationMode
{
    SynonymousOnly,
    NonsynonymousOnly,
    MutateAll
};

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;
static const int LATPACK_ERROR = 4;

// default values
static const int GENS_MAX = 99999;
static const std::string DEFAULT_OUTPATH = "./out";
static const uint64_t DEFAULT_SEED = 1;
static const int DEFAULT_CHECKPOINT_FREQ = 10;
static const int DEFAULT_JSON_OUTFREQ = 5;
static const int DEFAULT_POPULATION_SIZE = 500;
static const MutationMode DEFAULT_MUTATION_MODE = NonsynonymousOnly;
static const double DEFAULT_TEMPERATURE = 0.2;
// Not commandline options:
static const int DEFAULT_LATFOLD_OUTFREQ = 5000;
static const double DEFAULT_REEVALUATION_RATIO = 0.25;
static const double DEFAULT_FITNESS_CONSTANT = 0.25;
static const double DEFAULT_DEGRADATION_PARAM = 1000000;
static const double DEFAULT_SIM_TIME = 1e6;
static const double DEFAULT_CELL_TIME = 100e6;	   // Used to normalize fitness.

static const std::vector<Codon> STOP_CODONS = {N_UAA, N_UAG, N_UGA};

// For now these variables for MPI will be global
// until a better revision of the code
int g_world_rank, g_subcomm_rank, g_world_size, g_subcomm_size;
MPI_Status g_status;
MPI_Comm g_subcomm;


// --------------------------------------------------

// Function declarations

// Print this program's help message
void print_help();

// Print error message to rank 0
void print_error(const std::string& message, bool debug_mode);

// Check if given directory exists (From sodapop)
bool is_dir_exist(const std::string& path);

// Create directory from specified path (from sodapop)
bool make_path(const std::string& path);

// Remove path
// https://stackoverflow.com/questions/734717/how-to-delete-a-folder-in-c
void remove_dir(const char *path);

// Load simulation checkpoint file.
json open_checkpoint_file(const std::string& checkpoint_path);

// Save simulation state so that simulation can be resumed.
void write_checkpoint(
    const std::string& checkpoint_path,
    json& checkpoint);

// Convert a vector of strings into a vector of char pointers for use
// by execv family of commands. Adds a NULL element to the end.
std::vector<char *> string_vec_to_cstring_vec(
    std::vector<std::string>& string_vec);


// Build a vector of strings that are arguments to invoke
// latfold. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfold_command(
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    int sim_steps,
    double temperature,
    int output_frequency,
    double degradation_scale);


// Use vfork and exec to run a latfold simulation.
//
// @param latfold_command The vector of strings built
//        by compose_latfold_command.
// @param outfile latfold program output will go to this file.
void run_latfold(
    std::vector<std::string> latfold_command,
    const std::string & outfile);


// // Build a vector of strings that are arguments to invoke
// // latMapTraj. Note arguments -traj is not added here.
// std::vector<std::string> compose_latmaptraj_command(
//     const std::string& latpack_path,
//     const std::string& folded_conformation);


// // (this function is not currently used)
// // Use fork and exec to run n_simulations of latMapTraj to analyze the
// // latFoldVec simulations that were run.
// void run_latmaptraj(
//     std::vector<std::string> & latmaptraj_command,
//     const std::string& h5file_base,
//     int n_simulations=1);


// Analyze latfold simulations to average protein output.
// This is a parallel function.
//
// @param filename The partial name of the file hdf5 output was
//        written to.
// @param degradation_param Timescale that unfolded protein degradation
//        acts on.
// @param native_energy Value is updated with native energy.
double get_protein_output_avg(
    const std::string& filename,
    double* native_energy,
    double degradation_param,
    double* old_pnat_average,
    double* old_pnat_weight,
    double t_cell=DEFAULT_CELL_TIME);


// Calculate fitness from protein output.
//
// @param protein_output The fraction of successful folding
//        simulations
// @param f_0 Fitness function parameter
double calculate_fitness(
    double protein_output,
    double f_0=DEFAULT_FITNESS_CONSTANT);


// // Update protein output to incorporate new folding simulations.
// //
// // @param old_protein_output The previously estimated protein output.
// // @param new_protein_output The newly estimated protein output.
// // @param n_reevaluators The number of simulations run to estimate
// //        the new fitness value.
// // @param n_total The total number of CPUs running simulations.
// // @param n_gens_without_accept Number of generations since last
// //        accepted mutation.
// double reaverage_protein_output(
//     double old_protein_output,
//     double new_protein_output,
//     int n_reevaluators,
//     int n_total,
//     int n_gens_without_accept);


//
std::vector<int> mutate_sequence(
    const std::vector<int> & input_sequence,
    MutationMode mutation_mode);


// Print the header for simulation text log.
void print_header(
    std::ostream * outstream,
    const std::vector<AminoAcid> & aa_sequence,
    const std::vector<int> & nuc_sequence);


// Print out current state.
void print_state(
    std::ostream * outstream,
    int generation,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double native_energy,
    bool accept);


// Save current state to json log.
void save_state(
    json& json_log,
    int generation,
    int n_gens_without_accept,
    int mutation_type,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double old_protein_output,
    double new_protein_output,
    double old_pnat,
    double new_pnat,
    double native_energy,
    bool accept);


// Write json file.
void write_log(
    const json& json_log,
    const std::string& json_log_path);


// --------------------------------------------------
// --------------------== MAIN ==--------------------
// --------------------------------------------------
int main(int argc, char** argv)
{
    // MPI init stuff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_world_size);
    assert(g_world_size > 1);    // We need at least two CPUs

    // Split cpus into two groups
    // world_rank [0, reevaluation_size) does reevaluation
    int reevaluation_size = (int)g_world_size * DEFAULT_REEVALUATION_RATIO;
    bool proc_does_reevaluation = (g_world_rank < reevaluation_size);
    MPI_Comm_split(
	MPI_COMM_WORLD, proc_does_reevaluation, g_world_rank, &g_subcomm);
    MPI_Comm_rank(g_subcomm, &g_subcomm_rank);
    MPI_Comm_size(g_subcomm, &g_subcomm_size);
    assert(g_subcomm_size > 0);

    // Variables to be determined by commandline options: 

    // Required:
    std::string sequence;
    int n_gens = 0;			// How many generations to run

    // Other params:
    bool debug_mode = false;
    bool resume_from_checkpoint = false;
    int random_codons = 0;
    int population_size = DEFAULT_POPULATION_SIZE;
    std::string folded_conformation;
    std::string out_path = DEFAULT_OUTPATH;
    std::string checkpoint_path;
    std::string lat_sim_out_path;
    std::string json_log_path;	// we're going to log using json
    std::unordered_map<std::string, double> recorded_fitnesses;
    json checkpoint;
    json json_log;
    std::ostream * outstream = &std::cout;
    double degradation_param = DEFAULT_DEGRADATION_PARAM;
    double fitness_constant = DEFAULT_FITNESS_CONSTANT;
    uint64_t rng_seed = DEFAULT_SEED;
    MutationMode mutation_mode = DEFAULT_MUTATION_MODE;
    double temperature = DEFAULT_TEMPERATURE;
    // End variables to be determined by commandline options: 

    // These particular params currently don't have
    // commandline arguments
    int latfold_output_frequency = DEFAULT_LATFOLD_OUTFREQ;
    int checkpoint_frequency = DEFAULT_CHECKPOINT_FREQ;

    // Begin option handling (https://linux.die.net/man/3/getopt)
    int c, option_index = 0;
    /* 
     * A struct option has four fields:
     *  char *name: long name of option
     *  int has_arg: one of no_argument, required_argument, optional_argument
     *  int *flag, int val: 
     *   if flag is NULL, then val identifies the option
     *   else *flag will be set to val.
     *
     * Upon encountering an option, getopt_long returns val if flag is
     * NULL and 0 otherwise.
     */
    static struct option long_options[] = {
	{"activity-constant", required_argument, NULL, 'a'},
	{"from-checkpoint", required_argument, NULL, 'f'},
	{"degradation-param", required_argument, NULL, 'k'},
	{"mutation-mode", required_argument, NULL, 'm'},
	{"native-fold", required_argument, NULL, 'n'},
	{"new-ngens", required_argument, NULL, 2},
	{"random-codons", no_argument, &random_codons, 1},
	{"out-path", required_argument, NULL, 'o'},
	{"population-size", required_argument, NULL, 'p'},
	{"seed", required_argument, NULL, 'r'},
	{"temperature", required_argument, NULL, 't'},
	{"debug", no_argument, NULL, 'd'},
	{"help", no_argument, NULL, 'h'},
	{NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "a:f:k:m:n:o:p:r:t:dh",
			    long_options, &option_index))
	   != -1)
    {
	switch (c)
	{
	case 'h':
	    print_help();
	    exit(0);
	case 'f':
	    // We are resuming from checkpoint.
	    checkpoint_path = optarg;
	    resume_from_checkpoint = true;
	    break;
	case 'a':
	    try
	    {
		fitness_constant = std::stod(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert fitness constant: "
		    << optarg << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    debug_mode = true;
	    break;
	case 'd':
	    debug_mode = true;
	    break;
	case 0:
	    break;
	case 'k':
	    try
	    {
		degradation_param = std::stod(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert degradation-param to float: "
		    << optarg << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    break;
	case 'm':
	    try
	    {
		switch (std::stoi(optarg))
		{
		case 0:
		    // mutation_mode = SynonymousOnly;
		    mutation_mode = NonsynonymousOnly;
		    std::cerr << "Fitness saving means only mutation-mode 1 is supported"
			      << std::endl;
		    std::cerr << "Switching mutation mode to Nonsynonymousonly"
			      << std::endl;
		    break;
		case 1:
		    mutation_mode = NonsynonymousOnly;
		    break;
		case 2:
		    // mutation_mode = MutateAll;
		    mutation_mode = NonsynonymousOnly;
		    std::cerr << "Fitness saving means only mutation-mode 1 is supported"
			      << std::endl;
		    std::cerr << "Switching mutation mode to Nonsynonymousonly"
			      << std::endl;
		    break;
		default:
		    throw std::invalid_argument("Mutation mode not one of 0, 1, 2");
		}
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "--mutation-mode argument not one of 0, 1, 2: "
		    << optarg << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    break;
	case 'n':
	    folded_conformation = optarg;
	    break;
	case 2:
	    try
	    {
		n_gens = std::stoi(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert --new-ngens: " << optarg
		    << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }	    
	    break;
        case 'o':
            out_path = optarg;
            break;
	case 'p':
	    try
	    {
		population_size = std::stoi(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert --population-size: " << optarg
		    << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    break;
	case 'r':
	    try
	    {
		rng_seed = std::stoull(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert seed: " << optarg << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    break;
	case 't':
	    try
	    {
		temperature = std::stod(optarg);
	    }
	    catch (...)
	    {
		std::ostringstream err;
		err << "Failed to convert temp: " << optarg << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    break;
	case '?':
	    // getopt_long will print an error message
	    print_help();
	    exit(PARSE_ERROR);
	default:
	    // we shouldn't end up here i think
	    std::ostringstream err;
	    err << "getopt switch default; returned " << c << std::endl;
	    print_error(err.str(), debug_mode);
	    exit(PARSE_ERROR);
	}
    } // End while getopt_long

    // Continued parsing || checkpoint resume
    if (!resume_from_checkpoint)
    {
	// Positional arguments
	if (optind + 2 == argc)
	{
	    sequence = argv[optind++];
	    try
	    {
		n_gens = std::stoi(argv[optind]);
		if (n_gens > GENS_MAX)
		{
		    std::ostringstream err;
		    err << "Provided N_GENS argument " << n_gens
			<< " cannot exceed " << GENS_MAX << std::endl;
		    print_error(err.str(), debug_mode);
		    exit(PARSE_ERROR);
		}
	    }
	    catch (const std::invalid_argument& e)
	    {
		std::ostringstream err;
		err << "Provided N_GENS argument " << argv[optind]
		    << " could not be converted to uint" << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	    catch (const std::out_of_range& e)
	    {
		// probably unlikely
		std::ostringstream err;
		err << "Provided N_GENS argument " << argv[optind]
		    << " is out of range" << std::endl;
		print_error(err.str(), debug_mode);
		exit(PARSE_ERROR);
	    }
	}
	else
	{
	    std::ostringstream err;
	    err << "Expected two positional arguments: SEQUENCE NROUNDS"
		<< std::endl;
	    for (int i =0; i<argc; ++i)
		err << argv[i] << " ";
	    err << std::endl;
	    err << "optind value: " << optind << std::endl;
	    print_error(err.str(), debug_mode);
	    print_help();
	    exit(PARSE_ERROR);
	}

	// Complain about missing parameters here
	if (folded_conformation.empty())
	{
	    std::ostringstream err;
	    err << "Argument --native-fold (-n) is mandatory" << std::endl;
	    print_error(err.str(), debug_mode);
	    exit(PARSE_ERROR);
	}

	lat_sim_out_path = out_path + "/latfold_simulations";
	json_log_path = out_path + "/simulation_log.json";

	checkpoint["out path"] = out_path;
	checkpoint["lat sim path"] = lat_sim_out_path;
	checkpoint["json log path"] = json_log_path;
	checkpoint["sequence"] = sequence;

	// Setup output directories
	if (!g_world_rank && is_dir_exist(out_path))
	{
	    std::ostringstream err;
	    err << "Error, output path exists: " << out_path << std::endl;
	    print_error(err.str(), debug_mode);
	    exit(IO_ERROR);
	}

	if (!g_world_rank)
	{
	    if (!make_path(out_path))
	    {
		std::ostringstream err;
		err << "Error, could not make output path: " << out_path << std::endl;
		print_error(err.str(), debug_mode);
		MPI_Abort(MPI_COMM_WORLD, IO_ERROR);
	    }
	    if (!make_path(lat_sim_out_path))
	    {
		std::ostringstream err;
		err << "Error, could not make output path: " << lat_sim_out_path
		    << std::endl;
		print_error(err.str(), debug_mode);
		MPI_Abort(MPI_COMM_WORLD, IO_ERROR);
	    }
	}

	// Wait for everyone
	MPI_Barrier(MPI_COMM_WORLD);
	// Wait in case network file system is slow
	while (!is_dir_exist(lat_sim_out_path))
	{
	    sleep(1);
	}
    
	// Initialize RNG, giving different processes different seeds
	set_threefry_array(rng_seed, g_world_rank, proc_does_reevaluation, 0);
    }
    else			// checkpoint resume
    {
	int simulations_per_gen;
	int checkpoint_reevaluation_size;

	checkpoint = open_checkpoint_file(checkpoint_path);
	MPI_Barrier(MPI_COMM_WORLD);
	checkpoint.at("out path").get_to(out_path);
	checkpoint.at("json log path").get_to(json_log_path);
	checkpoint.at("lat sim path").get_to(lat_sim_out_path);
	// The original sequence
	checkpoint.at("sequence").get_to(sequence);

	// Access data stored in json log of simulation being resumed.
	// This restores the simulation trajectory log and lets us
	//   access a bunch of simulation parameters.
	std::ifstream json_log_file(json_log_path);
	json_log_file >> json_log;
	json_log.at("folded conformation").get_to(folded_conformation);
	if (n_gens == 0)
	{
	    json_log.at("generations").get_to(n_gens);
	}
	else
	{
	    json_log["generations"] = n_gens;
	}
	json_log.at("population size").get_to(population_size);
	json_log.at("temperature").get_to(temperature);
	json_log.at("simulations per gen").get_to(simulations_per_gen);
	json_log.at("reevaluation size").get_to(checkpoint_reevaluation_size);
	json_log.at("degradation timescale").get_to(degradation_param);
	json_log.at("fitness constant").get_to(fitness_constant);
	std::string latpack_share_path;
	std::string current_latpack_share_path(getenv("LATPACK_SHARE"));
	json_log.at("latpack share path").get_to(latpack_share_path);
	std::string foldevo_share_path;
	std::string current_foldevo_share_path(getenv("FOLDEVO_SHARE"));
	json_log.at("foldevo share path").get_to(foldevo_share_path);
	json_log.at("mutation mode").get_to(mutation_mode);

	if (latpack_share_path != current_latpack_share_path)
	{
	    std::cerr << "Current LATPACK_SHARE and checkpoint LATPACK_SHARE "
		      << "do not match.\n"
		      << current_latpack_share_path << std::endl
		      << " vs\n"
		      << latpack_share_path << std::endl;
	    std::cerr << "Use right latpack version.\n";
	    exit(DATA_ERROR);
	}

	if (foldevo_share_path != current_foldevo_share_path)
	{
	    std::cerr << "Current FOLDEVO_SHARE and checkpoint FOLDEVO_SHARE "
		      << "do not match.\n"
		      << current_foldevo_share_path << std::endl
		      << " vs\n"
		      << foldevo_share_path << std::endl;
	    std::cerr << "Use right foldevo version.\n";
	    exit(DATA_ERROR);
	}

	// Check that simulations_per_gen agrees with number of processors
	if (simulations_per_gen != g_world_size)
	{
	    std::cerr << "Resuming from simulation checkpoint, but number "
		      << "of CPUs currently does not match that of simulation "
		      << "that produced the checkpoint." << std::endl;
	    std::cerr << "Current CPU count: " << g_world_size << std::endl;
	    std::cerr << "Checkpoint CPU count: " << simulations_per_gen
		      << std::endl;
	    exit(DATA_ERROR);
	}
	if (checkpoint_reevaluation_size != reevaluation_size)
	{
	    std::cerr << "Resuming from simulation checkpoint, but number "
		      << "of CPUs to be used to reevaluate accepted "
		      << "sequences does not match that of simulation "
		      << "that produced the checkpoint." << std::endl;
	    std::cerr << "Current reevaluation count: " << reevaluation_size
		      << std::endl;
	    std::cerr << "Checkpoint reevaluation count: "
		      << checkpoint_reevaluation_size << std::endl;
	    exit(DATA_ERROR);
	}

	// need to trim json log if it's too long
	// also, we will reconstruct recorded_fitnesses
	int checkpoint_generation;
	checkpoint.at("generation").get_to(checkpoint_generation);
	auto it = json_log.at("trajectory").begin();
	std::string curr_seq = it->at("aa sequence");
	std::string new_seq;
	while (it != json_log.at("trajectory").end() &&
	       it->at("generation") <= checkpoint_generation)
	{
	    new_seq = it->at("aa sequence");
	    if (it->at("accepted"))
	    {
		recorded_fitnesses[curr_seq] = it->at("old fitness");
		curr_seq = new_seq;
	    }
	    else
	    {
		recorded_fitnesses[new_seq] = it->at("new fitness");
	    }
	    ++it;
	}

	// auto it2 = it;
	// for (; it2 != json_log.at("trajectory").end(); ++it2)
	// {
	//     recorded_fitnesses.erase(it2->at("aa sequence"));
	// }

	json_log.at("trajectory").erase(it, json_log.at("trajectory").end());

    }
    
    // End option parsing
    // --------------------------------------------------
    // Begin simulation setup

    // Process the user-provided sequence
    std::vector<int> nuc_sequence;
    std::vector<int> old_nuc_sequence;
    std::vector<AminoAcid> aa_sequence;
    std::vector<AminoAcid> old_aa_sequence;
    std::vector<Codon> codon_sequence;

    // Infer lengths from folded conformation length
    unsigned int protein_length = folded_conformation.length() + 1;
    unsigned int nuc_length = protein_length * 3;

    // Resize vectors
    nuc_sequence.resize(nuc_length);
    old_nuc_sequence.resize(nuc_length);
    aa_sequence.resize(protein_length);
    old_aa_sequence.resize(protein_length); // for recorded_fitnesses
    codon_sequence.resize(protein_length);

    if (!resume_from_checkpoint)
    {
	// We use resize because we directly access the underlying arrays of
	// our three vectors in this next section (b/c interfacing with C code)
	if (protein_length == sequence.length()) 
	{
	    // Presume it's a protein sequence
	    if (LetterToAASeq(
		    sequence.c_str(), aa_sequence.data(), protein_length))
	    {
		std::ostringstream err;
		err << "Letter to aa conversion problem" << std::endl;
		err << "Sequence: " << sequence << std::endl;
		print_error(err.str(), debug_mode);
		exit(DATA_ERROR);
	    }
	    AASeqToNucSeq(aa_sequence.data(), nuc_sequence.data(),
			  protein_length, random_codons);
	}
	else if (protein_length * 3 == sequence.length())
	{
	    // presume we have an rna sequence
	    if (LetterToNucCodeSeq(
		    sequence.c_str(), nuc_sequence.data(), nuc_length))
	    {
		std::ostringstream err;
		err << "letter to nuc conversion problem" << std::endl;
		err << "Sequence: " << sequence << std::endl;
		print_error(err.str(), debug_mode);
		exit(DATA_ERROR);
	    }
	    if (NucSeqToAASeq(
		    nuc_sequence.data(), nuc_length, aa_sequence.data()))
	    {
		// nuc sequence contains a stop codon
		std::ostringstream err;
		err << "Provided RNA sequence contains a stop codon"
		    << std::endl;
		err << "Sequence: " << sequence << std::endl;
		print_error(err.str(), debug_mode);
		exit(DATA_ERROR);
	    }
	}
	else
	{
	    // houston we have a problem
	    std::ostringstream err;
	    err << "Protein length of " << protein_length
		<< ", inferred from conformation string," << std::endl
		<< "is not matched by provided sequence:" << std::endl
		<< sequence << ", l = " << sequence.length() << std::endl
		<< "Expected length of " << protein_length << " or "
		<< nuc_length << " (RNA sequence)" << std::endl;
	    print_error(err.str(), debug_mode);
	    exit(DATA_ERROR);
	}
    }
    else			// resuming from checkpoint
    {
	nuc_sequence.clear();
	old_nuc_sequence.clear();
	checkpoint.at("nuc sequence").get_to(nuc_sequence);
	checkpoint.at("old nuc sequence").get_to(old_nuc_sequence);
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());
    }

    // Obtain a codon sequence from nuc_sequence
    NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());

    // It's time to give the people some information
    if (!g_world_rank)
    {
	if (!resume_from_checkpoint)
	{
	    json_log["sequence"] = sequence;
	    json_log["folded conformation"] = folded_conformation;
	    json_log["generations"] = n_gens;
	    json_log["population size"] = population_size;
	    json_log["temperature"] = temperature;
	    json_log["simulations per gen"] = g_world_size;
	    json_log["reevaluation size"] = reevaluation_size;
	    json_log["degradation timescale"] = degradation_param;
	    json_log["fitness constant"] = fitness_constant;
	    json_log["latpack share path"] = getenv("LATPACK_SHARE");
	    json_log["foldevo share path"] = getenv("FOLDEVO_SHARE");
	    json_log["rng seed"] = rng_seed;
	    json_log["mutation mode"] = mutation_mode;
	    json_log["trajectory"] = json::array();
	}
	// printing
	*outstream
	    << "# folding_evolution" << std::endl
	    << "# input seq : " << sequence << std::endl
	    << "# n_gens : " << n_gens << std::endl
	    << "# pop size : " << population_size << std::endl
	    << "# temperature : " << temperature << std::endl
	    << "# simulations (no. reevaluations): " << g_world_size <<
	    " (" << reevaluation_size << ")" << std::endl
	    << "# degradation timescale : " << degradation_param << std::endl
	    << "# fitness constant : " << fitness_constant << std::endl
	    << "# rng seed : " << rng_seed << std::endl
	    << "# latpack share : " << getenv("LATPACK_SHARE") << std::endl
	    ;
	print_header(outstream, aa_sequence, nuc_sequence); 
    }
    
    int folding_time = DEFAULT_SIM_TIME * protein_length;
    double degradation_scale_steps = degradation_param * protein_length;
    int gen = 0;
    int last_accepted_gen = 0;
    int n_gens_without_accept = 0;
    int mutation_type = -1;
    double old_fitness = 0.0000001;
    double fitness;
    double old_protein_output = 0;
    double protein_output;
    double pnat_average;
    double old_pnat_average = 0;
    double pnat_weight;
    double old_pnat_weight;
    double native_energy;
    double selection;
    double fixation;
    std::vector<std::string> latfold_command;
    std::vector<std::string> old_latfold_command;
    std::unique_ptr<char []> aa_sequence_str(new char[protein_length+1]); // the current trial sequence
    std::unique_ptr<char []> old_aa_sequence_str(new char[protein_length+1]); // the current accepted

    if (resume_from_checkpoint)
    {
	checkpoint.at("generation").get_to(gen);
	++gen;
	checkpoint.at("last accepted gen").get_to(last_accepted_gen);
	checkpoint.at("n gens without accept").get_to(
	    n_gens_without_accept);
	checkpoint.at("mutation type").get_to(mutation_type);
	checkpoint.at("old fitness").get_to(old_fitness);
	checkpoint.at("old protein output").get_to(old_protein_output);
	checkpoint.at("old latfold command").get_to(
	    old_latfold_command);
	checkpoint.at("old pnat average").get_to(old_pnat_average);
	checkpoint.at("old pnat weight").get_to(old_pnat_weight);
    }

    // Begin running simulation loop
    for (; gen <= n_gens; ++gen)
    {
	PrintAASequence(
	    aa_sequence_str.get(), aa_sequence.data(), protein_length);
	latfold_command = compose_latfold_command(
	    aa_sequence_str.get(),
	    folded_conformation,
	    folding_time,
	    temperature,
	    latfold_output_frequency,
	    degradation_scale_steps);

	// Reset these variables for new eval
	pnat_average = 0;
	pnat_weight = 0;

	if (proc_does_reevaluation && gen > 0)
	{
	    std::ostringstream output_dir;
	    std::ostringstream hdf5_output_file;
	    output_dir << lat_sim_out_path << "/gen" << std::setw(5)
		       << std::setfill('0') << last_accepted_gen;
	    hdf5_output_file << output_dir.str() << "/sim" << std::setfill('0')
			     << std::setw(5) << n_gens_without_accept + 1
			     << "_" << std::setw(5) << g_subcomm_rank << ".h5";

	    run_latfold(old_latfold_command, hdf5_output_file.str());
	    old_protein_output = get_protein_output_avg(
		hdf5_output_file.str(), &native_energy,
		degradation_param, &old_pnat_average, &old_pnat_weight);
	    old_fitness = calculate_fitness(
		old_protein_output, fitness_constant);
	}
	else if (!proc_does_reevaluation)
	{
	    std::ostringstream output_dir;
	    std::ostringstream hdf5_output_file;
	    output_dir << lat_sim_out_path << "/gen" << std::setw(5)
		       << std::setfill('0') << std::to_string(gen);

	    if (!g_subcomm_rank)
	    {
		if (!make_path(output_dir.str()))
		{
		    std::cerr << "Error making dir: " << output_dir.str()
			      << std::endl;
		    MPI_Abort(MPI_COMM_WORLD, IO_ERROR);
		}
	    }

            // wait for mkdir
	    MPI_Barrier(g_subcomm); 
	    while (!is_dir_exist(output_dir.str()))
	    {
		sleep(1);
	    }

	    hdf5_output_file << output_dir.str() << "/sim" << std::setfill('0')
			     << std::setw(5) << 0
			     << "_" << std::setw(5) << g_subcomm_rank << ".h5";

	    run_latfold(latfold_command, hdf5_output_file.str());
	    protein_output = get_protein_output_avg(
		hdf5_output_file.str(), &native_energy,
		degradation_param, &pnat_average, &pnat_weight);
	    fitness = calculate_fitness(
		protein_output, fitness_constant);
	}

	// Synchronize data between all ranks, including between
	// reevaluators and evaluators. ranks [0, reevaluation_size)
	// looks at old sequence whereas [reevaluation_size,
	// g_world_size) does new evaluation
	MPI_Bcast(&native_energy, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&fitness, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&protein_output, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&pnat_average, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&pnat_weight, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&old_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&old_protein_output, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&old_pnat_average, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&old_pnat_weight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	selection = (fitness / old_fitness) - 1;
	if (selection == 0.)
	    fixation = 1 / (double)population_size;
	else
	    fixation = (1 - exp(-2 * selection)) /
		(1 - exp(-2 * population_size * selection));

	// Decide whether to accept
	int accept;
	if (!g_world_rank)
	{
	    // Root node decide accept / reject
	    if (threefryrand() < fixation || gen == 0)
		accept = true;
	    else
		accept = false;

	    // Output information
	    print_state(outstream, gen, aa_sequence, nuc_sequence, old_fitness,
			fitness, native_energy, accept);
	    save_state(json_log, gen, n_gens_without_accept, mutation_type,
		       aa_sequence, nuc_sequence, old_fitness, fitness,
		       old_protein_output, protein_output, old_pnat_average,
		       pnat_average, native_energy, (bool)accept);
	    if (gen % DEFAULT_JSON_OUTFREQ == 0)
	    {
		write_log(json_log, json_log_path);
	    }
	}

	// Let everyone know accept / reject
	MPI_Bcast(&accept, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	// Update / revert
	if (accept)
	{
	    if (gen > 0)
	    {
		NucSeqToAASeq(old_nuc_sequence.data(), nuc_length, old_aa_sequence.data());
		PrintAASequence(
		    old_aa_sequence_str.get(), old_aa_sequence.data(), protein_length);
		recorded_fitnesses[old_aa_sequence_str.get()] = old_fitness;
	    }
	    // Update old values to take on current sequence
	    old_nuc_sequence = nuc_sequence;
	    old_latfold_command = latfold_command;
	    old_protein_output = protein_output;
	    old_fitness = fitness;
	    n_gens_without_accept = 0;
	    last_accepted_gen = gen;
	    old_pnat_average = pnat_average;
	    old_pnat_weight = pnat_weight;
	}
	else
	{
	    // Revert, saving the fitness value
	    recorded_fitnesses[aa_sequence_str.get()] = fitness;
	    nuc_sequence = old_nuc_sequence;
	    n_gens_without_accept++;
	}

	// Make new mutation. We only do this on root node because of
	// independent rng streams for different nodes.
	if (!g_world_rank)
	{
	    // We check if mutant sequence has been encountered before.
	    std::vector<int> new_nuc_sequence;
	    std::vector<AminoAcid> new_aa_sequence;
	    new_aa_sequence.resize(protein_length);
	    int infinite_loop_checker = 0;
	    while (true)
	    {
		new_nuc_sequence = mutate_sequence(nuc_sequence, mutation_mode);
		NucSeqToAASeq(new_nuc_sequence.data(), nuc_length, new_aa_sequence.data());
		PrintAASequence(
		    aa_sequence_str.get(), new_aa_sequence.data(), protein_length);
		auto result = recorded_fitnesses.find(aa_sequence_str.get());
		if (result == recorded_fitnesses.end())
		{
		    break;
		}
		else
		{
		    if (result->second >= old_fitness)
		    {
			break;
		    }
		}

		++infinite_loop_checker;
		if (infinite_loop_checker > 10000)
		{
		    std::cerr << "# No new sequences / s>=0 sequences in 10000 mutation attempts. "
			      << "Fitness maxima? Continuing with current mutant.\n";
		    break;
		}
	    }
	    nuc_sequence = new_nuc_sequence;
	}

	MPI_Bcast(nuc_sequence.data(), nuc_sequence.size(), MPI_INT,
		  0, MPI_COMM_WORLD);
	
	// and update sequences
	NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());

	// Need to save checkpoint
	if (gen % checkpoint_frequency == 0 || gen == n_gens)
	{
	    std::ostringstream checkpoint_path;
	    checkpoint_path << out_path << "/checkpoint"<< std::setw(5)
			   << std::setfill('0') << std::to_string(gen)
			   << ".json";

	    checkpoint["generation"] = gen;
	    checkpoint["nuc sequence"] = nuc_sequence;
	    checkpoint["old nuc sequence"] = old_nuc_sequence;
	    checkpoint["last accepted gen"] = last_accepted_gen;
	    checkpoint["n gens without accept"] = n_gens_without_accept;
	    checkpoint["mutation type"] = mutation_type;
	    checkpoint["old fitness"] = old_fitness;
	    checkpoint["old protein output"] = old_protein_output;
	    checkpoint["old latfold command"] = old_latfold_command;
	    checkpoint["old pnat average"] = old_pnat_average;
	    checkpoint["old pnat weight"] = old_pnat_weight;
	    checkpoint["recorded fitnesses"] = recorded_fitnesses;
	    write_checkpoint(checkpoint_path.str(), checkpoint);
	}
    }

    // finalize json log
    write_log(json_log, json_log_path);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Comm_free(&g_subcomm);
    MPI_Finalize();

    return 0;
} // end main


// Print this program's help message
void print_help()
{
    if (!g_world_rank)
	std::cerr << helptext;
}


// Print error message to rank 0
void print_error(const std::string& message, bool debug_mode)
{
    if (!g_world_rank || debug_mode)
	std::cerr << message;
}

// checks if given directory exists
bool is_dir_exist(const std::string& path)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
        return false;
    return (info.st_mode & S_IFDIR) != 0;
}


// creates directory from specified path
bool make_path(const std::string& path)
{
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
    if (ret == 0)
        return true;

    // else...
    switch (errno)
    {
    case ENOENT:
        // parent didn't exist, try to create it
    {
	std::size_t pos = path.find_last_of('/');
	if (pos == std::string::npos)
	    return false;
	if (!make_path( path.substr(0, pos)))
	    return false;
    }
    return 0 == mkdir(path.c_str(), mode);
    ;
    case EEXIST:
        // done
        return is_dir_exist(path);

    default:
        return false;
    }
}


// https://stackoverflow.com/questions/734717/how-to-delete-a-folder-in-c
void remove_dir(const char *path)
{
    struct dirent *entry = NULL;
    DIR *dir = NULL;
    dir = opendir(path);
    while ((entry = readdir(dir)))
    {   
	DIR *sub_dir = NULL;
	FILE *file = NULL;
	char* abs_path = new char[1024];
	if ((*(entry->d_name) != '.') || ((strlen(entry->d_name) > 1) && (entry->d_name[1] != '.')))
	{   
	    sprintf(abs_path, "%s/%s", path, entry->d_name);
	    if((sub_dir = opendir(abs_path)))
	    {   
		closedir(sub_dir);
		remove_dir(abs_path);
	    }   
	    else 
	    {   
		if((file = fopen(abs_path, "r")))
		{   
		    fclose(file);
		    remove(abs_path);
		}   
	    }   
	}
	delete[] abs_path;   
    }   
    remove(path);
}


// Load simulation checkpoint file.
json open_checkpoint_file(const std::string & checkpoint_path)
{
    std::ifstream infile(checkpoint_path);
    json checkpoint;
    infile >> checkpoint;

    // Load RNG state.
    json rng_keys = checkpoint.at("rng key");
    json rng_counter_values = checkpoint.at("rng counter");
    json rng_result_values = checkpoint.at("rng result");
    json rng_indices = checkpoint.at("rng index");
    json counter_indices = checkpoint.at("rng counter index");
    std::vector<uint64_t> rng_key;
    std::vector<uint64_t> rng_counter;
    std::vector<uint64_t> rng_result;
    unsigned char rng_index;
    unsigned char counter_index;
    rng_keys[g_world_rank].get_to(rng_key);
    rng_counter_values[g_world_rank].get_to(rng_counter);
    rng_result_values[g_world_rank].get_to(rng_result);
    rng_indices[g_world_rank].get_to(rng_index);
    counter_indices[g_world_rank].get_to(counter_index);
    set_threefry_array(rng_key[0], rng_key[1], rng_key[2], rng_key[3]);
    set_threefry_counter(rng_counter[0], rng_counter[1],
			 rng_counter[2], rng_counter[3], counter_index);
    set_threefry_result(rng_result[0], rng_result[1],
			rng_result[2], rng_result[3], rng_index);


    // Delete leftover files:
    if (!g_world_rank)
    {
	int checkpoint_gen;
	std::string lat_sim_out_path;
	checkpoint.at("generation").get_to(checkpoint_gen);
	checkpoint.at("lat sim path").get_to(lat_sim_out_path);
    
	// Delete generations after checkpoint:
	std::string pathway_base(lat_sim_out_path + "/gen");
	std::ostringstream pathname;
	int gen = checkpoint_gen + 1;
	while (true)
	{
	    pathname << lat_sim_out_path + "/gen" << std::setw(5)
		     << std::setfill('0') << std::to_string(gen);
	    if (access(pathname.str().c_str(), F_OK) != -1)
	    {
		// Remove files.
		remove_dir(pathname.str().c_str());
	    }
	    else
	    {
		break;
	    }
	    pathname.str("");
	    gen++;
	}
    }

    return checkpoint;
}


// Save simulation state so that simulation can be resumed.
void write_checkpoint(
    const std::string& checkpoint_path,
    json& checkpoint)
{
    // things that need to be saved:
    // generation number, random number of each processor, number of processors,

    // Save RNG state
    checkpoint["rng key"] = json::array();
    checkpoint["rng counter"] = json::array();
    checkpoint["rng result"] = json::array();
    checkpoint["rng index"] = json::array();
    checkpoint["rng counter index"] = json::array();
    std::vector<uint64_t> rng_key(4, 0);
    std::vector<uint64_t> rng_counter(4, 0);
    std::vector<uint64_t> rng_result(4, 0);
    unsigned char rng_index;
    unsigned char counter_index;
    get_rng_state(rng_key.data(), rng_counter.data(), rng_result.data(),
		  &rng_index, &counter_index);
    if (g_world_rank)
    {
	// mpi send
	MPI_Send(rng_key.data(), 4, MPI_UINT64_T, 0, 0, MPI_COMM_WORLD);
	MPI_Send(rng_counter.data(), 4, MPI_UINT64_T, 0, 1, MPI_COMM_WORLD);
	MPI_Send(rng_result.data(), 4, MPI_UINT64_T, 0, 2, MPI_COMM_WORLD);
	MPI_Send(&rng_index, 1, MPI_UINT8_T, 0, 3, MPI_COMM_WORLD);
	MPI_Send(&counter_index, 1, MPI_UINT8_T, 0, 4, MPI_COMM_WORLD);
    }
    else
    {
	checkpoint["rng key"].push_back(rng_key);
	checkpoint["rng counter"].push_back(rng_counter);
	checkpoint["rng result"].push_back(rng_result);
	checkpoint["rng index"].push_back(rng_index);
	checkpoint["rng counter index"].push_back(counter_index);
	// mpi recv
	for (int i=1; i<g_world_size; ++i)
	{
	    MPI_Recv(rng_key.data(), 4, MPI_UINT64_T, i, 0, MPI_COMM_WORLD,
		&g_status);
	    MPI_Recv(rng_counter.data(), 4, MPI_UINT64_T, i, 1, MPI_COMM_WORLD,
		&g_status);
	    MPI_Recv(rng_result.data(), 4, MPI_UINT64_T, i, 2, MPI_COMM_WORLD,
		&g_status);
	    MPI_Recv(&rng_index, 1, MPI_UINT8_T, i, 3, MPI_COMM_WORLD,
		&g_status);
	    MPI_Recv(&counter_index, 1, MPI_UINT8_T, i, 4, MPI_COMM_WORLD,
		&g_status);	    
	    checkpoint["rng key"].push_back(rng_key);
	    checkpoint["rng counter"].push_back(rng_counter);
	    checkpoint["rng result"].push_back(rng_result);
	    checkpoint["rng index"].push_back(rng_index);
	    checkpoint["rng counter index"].push_back(counter_index);
	}
    }

    // Save from root.
    if (!g_world_rank)
    {
	std::ofstream checkpointFile(checkpoint_path);
	checkpointFile << std::setw(2) << checkpoint << std::endl;
    }
    return;
}


// Convert a vector of strings into a vector of char pointers for use
// by execv family of commands. Adds a NULL element to the end.
std::vector<char *> string_vec_to_cstring_vec(
    std::vector<std::string>& string_vec)
{
    std::vector<char *> cstrings;
    for (auto& string : string_vec)
    {
	cstrings.push_back(&string.front());
    }
    cstrings.push_back(NULL);
    return cstrings;
}


// Build a vector of strings that are arguments to invoke
// latfold. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfold_command(
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    int sim_steps,
    double temperature,
    int output_frequency,
    double degradation_scale)	// in steps
{
    std::vector<std::string> command;

    command.push_back("latFold");
    std::string energyFileArg(getenv("LATPACK_SHARE"));
    command.push_back("-energyFile=" + energyFileArg + "/MJ.txt");
    command.push_back("-seq=" + aa_sequence);
    command.push_back("-abs=" + folded_conformation);
    command.push_back("-countTarget=" + folded_conformation);
    command.push_back("-targetDegradationScale=" + std::to_string(degradation_scale));
    command.push_back("-kT=" + std::to_string(temperature));
    command.push_back("-maxSteps=" + std::to_string(sim_steps));
    command.push_back("-outFreq=" + std::to_string(output_frequency));
    command.push_back("-out=S");
    // if (save_conformations)
    // {
    //     command.push_back("-out=S");
    // }
    // else
    // {
    // 	command.push_back("-out=E");
    // }
    
    // Previously we did the what is commented out, but we want to
    // keep the command as a vector.
    // std::ostringstream os;
    // std::copy(command.begin(), command.end(),
    // std::ostream_iterator<std::string>(os, " ")); return os.str();
    return command;
}


// Use fork and exec to run a latfold co-translational folding
// simulation.
void run_latfold(
    std::vector<std::string> latfold_command,
    const std::string & outfile)
{
    // Add additional parameters
    latfold_command.push_back(
	"-seed=" + std::to_string(threefryrand_int() % INT32_MAX));
    latfold_command.push_back("-outFile=" + outfile);
    // latfold_command.push_back("-title=?");

    std::ostringstream command;
    for (auto& arg : latfold_command)
    {
	command << arg << " ";
    }

    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
	latfold_command);

    // Temporarily suppress stdout
    int bak, temp_fd;
    fflush(stdout);
    bak = dup(1);		// backup
    temp_fd = open("/dev/null", O_WRONLY);
    dup2(temp_fd, 1);
    close(temp_fd);

    pid_t pid = vfork();
    if (pid == -1)
    {
	std::cerr << "Failed to fork" << std::endl;
	exit(IO_ERROR);
    }
    else if (pid == 0)
    {
	// child
	// suppress std out
	// int fd = open("/dev/null", O_WRONLY);
	// dup2(fd, 1);

	execvp(cstring_command_vec[0], cstring_command_vec.data());
    }

    // wait for processes to finish
    while (true)
    {
	int status;
	pid_t done = wait(&status);
	if (done == -1)
	{
	    if (errno == ECHILD) break; // no more child processes
	}
	else
	{
	    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
	    {
		std::cerr << "pid " << done << " failed" << std::endl;
		std::cerr << command.str() << std::endl;
		exit(LATPACK_ERROR);
	    }
	}
    }

    // Restore stdout
    fflush(stdout);
    dup2(bak, 1);
    close(bak);
    return;
}


// Build a vector of strings that are arguments to invoke
// latMapTraj. Note arguments -traj is not added here.
// std::vector<std::string> compose_latmaptraj_command(
//     const std::string& latpack_path,
//     const std::string& folded_conformation)
// {
//     std::vector<std::string> command;
//     command.push_back(latpack_path + "bin/latMapTraj");
//     command.push_back("-ref=" + folded_conformation);
//     command.push_back("-lat=CUB");
//     command.push_back("-hdf5");
//     return command;
// }


// FIXME (update to full MPI)
// Use fork and exec to run n_simulations of latMapTraj to analyze the
// latfold simulations that were run.
// this has not been updated yet (because function is not used)
// void run_latmaptraj(
//     std::vector<std::string> & latmaptraj_command,
//     const std::string& h5file_base,
//     int n_simulations)
// {
//     for (int i=g_world_rank*n_simulations; i<(g_world_rank+1)*n_simulations; ++i)
//     {
// 	pid_t pid = fork();
// 	if (pid == -1)
// 	{
// 	    std::cerr << "Failed to fork" << std::endl;
// 	    exit(IO_ERROR);
// 	}
// 	else if (pid == 0)
// 	{
// 	    // child
// 	    // we land here for each instance
// 	    // Add additional parameters
// 	    latmaptraj_command.push_back("-traj=" + h5file_base +
// 					 std::to_string(i) + ".h5");
// 	    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
// 		latmaptraj_command);

// 	    // suppress std out
// 	    int fd = open("/dev/null", O_WRONLY);
// 	    dup2(fd, 1);

// 	    execv(cstring_command_vec[0], cstring_command_vec.data());
// 	}
//     }

//     // wait for processes to finish
//     while (true)
//     {
// 	int status;
// 	pid_t done = wait(&status);
// 	if (done == -1)
// 	{
// 	    if (errno == ECHILD) break; // no more child processes
// 	}
// 	else
// 	{
// 	    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
// 	    {
// 		std::cerr << "pid " << done << " failed" << std::endl;
// 		exit(LATPACK_ERROR);
// 	    }
// 	}
//     }    
// }



// Analyze latfold simulations to average protein output.
// This is a parallel function.
// This version of code uses pnat only
double get_protein_output_avg(
    const std::string& filename,
    double* native_energy,
    double degradation_param,
    double* old_pnat_average,
    double* old_pnat_weight,
    double t_cell)
{
    // We need to read data from HDF5 file
    int error_count = 0;
    hid_t file_id;
    while (error_count < 5)
    {
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id >= 0)
	    break;
	sleep(5);
	++error_count;
    }
    if (error_count > 5)
    {
	std::cerr << "Could not open file \"" << filename
		  << "\" after 5 attempts.\n";
	std::cerr << "Ending program.\n";
	MPI_Abort(MPI_COMM_WORLD, IO_ERROR);	
    }

    hid_t group_id = H5Gopen2(file_id, "traj1", H5P_DEFAULT);

    // Read attributes from 'traj1'
    hid_t pnat_attr = H5Aopen(group_id, "Fraction target conf",
			      H5P_DEFAULT);
    hid_t conf_energy_attr = H5Aopen(group_id, "Target conf energy",
				     H5P_DEFAULT);
    double pnat = 0;
    double conf_energy = 0;

    H5Aread(pnat_attr, H5T_NATIVE_DOUBLE, &pnat);
    H5Aread(conf_energy_attr, H5T_NATIVE_DOUBLE, &conf_energy);

    H5Aclose(conf_energy_attr);
    H5Aclose(pnat_attr);
    H5Gclose(group_id);
    H5Fclose(file_id);
    // End reading hdf5 file

    // Get a weighted average for pnat. Gather everything on root.
    std::vector<double> pnats(g_subcomm_size);

    MPI_Gather(&pnat, 1, MPI_DOUBLE, pnats.data(),
	       1, MPI_DOUBLE, 0, g_subcomm);

    double output_average;
    if (g_subcomm_rank == 0)
    {
	double sum = 0;
	double sum_weights = 0;
	for (unsigned int i = 0; i < pnats.size(); ++i)
	{
	    sum += pnats[i];
	    ++sum_weights;
	}
	pnat = sum / sum_weights;

	if (std::isnan(pnat))
	{
	    pnat = 0;
	}

	// Now account for old:
	pnat = pnat * sum_weights + *old_pnat_average * *old_pnat_weight;
	*old_pnat_weight += sum_weights;
	pnat /= *old_pnat_weight;
	*old_pnat_average = pnat;

	if (pnat == 1)
	{
	    pnat = 1 - 1e-10;
	}

	output_average = degradation_param * pnat / (1 - pnat);
	output_average *= (1 - exp(-t_cell * (1 - pnat)
				   / degradation_param)) / t_cell;
    }

    // Now let everyone know
    MPI_Bcast(old_pnat_average, 1, MPI_DOUBLE, 0, g_subcomm);
    MPI_Bcast(old_pnat_weight, 1, MPI_DOUBLE, 0, g_subcomm);
    MPI_Bcast(&output_average, 1, MPI_DOUBLE, 0, g_subcomm);
    pnat = *old_pnat_average;

    // Also determine native energy
    MPI_Allreduce(&conf_energy, native_energy, 1, MPI_DOUBLE, MPI_MIN,
		  g_subcomm);
    
    return output_average;
}


// Calculate fitness from protein output.
double calculate_fitness(
    double protein_output,
    double f_0)
{
    // 0.0000001 (1e-7) to avoid divide by zero error
    return 0.0000001 + protein_output / (protein_output + f_0);
    // return 0.0000001 + protein_output;
}


// // Update protein output to incorporate new folding simulations.
// double reaverage_protein_output(
//     double old_protein_output,
//     double new_protein_output,
//     int n_reevaluators,
//     int n_total,
//     int n_gens_without_accept)
// {
//     // We calculate a weighted average, so here we get the weights.
//     double old_factor = (n_total + n_reevaluators * (n_gens_without_accept - 1))
// 	/ (double)(n_total + n_gens_without_accept * n_reevaluators);
//     double new_factor = n_reevaluators
// 	/ (double)(n_total + n_gens_without_accept * n_reevaluators);
//     return old_protein_output * old_factor + new_protein_output * new_factor;
// }


// Make a mutation
std::vector<int> mutate_sequence(
    const std::vector<int> & input_sequence,
    MutationMode mutation_mode)
{
    std::vector<int> temp_sequence;
    int mutation_type;
    unsigned int nuc_length = input_sequence.size();
    if (mutation_mode == NonsynonymousOnly)
    {
	temp_sequence = input_sequence;
	AAMutateNucSequence(temp_sequence.data(), nuc_length);
	// This function should always work.
    }
    else		// synonymousonly == 0, mutateall == 2
    {
	while (true)
	{
	    temp_sequence = input_sequence;
	    mutation_type = PointMutateNucSequence(
		temp_sequence.data(), nuc_length);
	    if (mutation_type == mutation_mode)
		// mutation_type == mutation_mode == 0
		break;
	    if (mutation_mode == MutateAll && mutation_type >= 0)
		// either
		// mutation_type == 0
		// mutation_type == 1
		break;
	}
    }
    return temp_sequence;
}


// Print the header for state output
void print_header(
    std::ostream * outstream,
    const std::vector<AminoAcid> & aa_sequence,
    const std::vector<int> & nuc_sequence)
{
    auto shorten = [outstream](std::string str, int len)
    {
	*outstream << std::setw(len);
	return str.substr(0, len);
    };

    *outstream << std::endl;
    *outstream << std::left <<	"# Gen|";
    *outstream << shorten("AA sequence", aa_sequence.size()) << "|";
    *outstream << shorten("Nuc sequence", nuc_sequence.size()) << "|";
    *outstream << 
	std::setw(8) << "Old fit." << "|" <<
	std::setw(8) << "New fit." << "|" <<
	std::setw(7) << "Nat. E." << "|" <<
	std::setw(6) << "Accept"
	       << std::endl;

    int line_length = 4;
    line_length += aa_sequence.size() + 1;
    line_length += nuc_sequence.size() + 1;
    line_length += 23 + 3;
    line_length += 6;

    *outstream << "# " << std::string(line_length, '-')<< std::endl;

    *outstream << std::fixed << std::setprecision(3);
    
    return;
}


// Print out what happened.
void print_state(
    std::ostream * outstream,
    int generation,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double native_energy,
    bool accept)
{
    auto aa_seq_len = aa_sequence.size();
    auto nuc_seq_len = nuc_sequence.size();
    std::unique_ptr<char []> aa_seq_str(new char[aa_seq_len+1]);
    std::unique_ptr<char []> nuc_seq_str(new char[nuc_seq_len+1]);
    PrintAASequence(aa_seq_str.get(), aa_sequence.data(), aa_seq_len);
    PrintNucCodeSequence(nuc_seq_str.get(), nuc_sequence.data(), nuc_seq_len);

    *outstream << std::right
	       << std::setw(5) << generation << " "
	       << std::setw(aa_seq_len) << aa_seq_str.get() << " "
	       << std::setw(nuc_seq_len) << nuc_seq_str.get() << " "
	       << std::setprecision(5)
	       << std::setw(8) << old_fitness << " "
	       << std::setw(8) << new_fitness << " "
	       << std::setprecision(2)
	       << std::setw(7) << native_energy << " "	
	       << std::setw(6) << ((accept) ? "yes" : "no")
	       << std::endl;

    return;    
}


// Save results to json log.
void save_state(
    json& json_log,
    int generation,
    int n_gens_without_accept,
    int mutation_type,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double old_protein_output,
    double new_protein_output,
    double old_pnat,
    double new_pnat,
    double native_energy,
    bool accept)
{
    auto aa_seq_len = aa_sequence.size();
    auto nuc_seq_len = nuc_sequence.size();
    std::unique_ptr<char []> aa_seq_str(new char[aa_seq_len+1]);
    std::unique_ptr<char []> nuc_seq_str(new char[nuc_seq_len+1]);
    PrintAASequence(aa_seq_str.get(), aa_sequence.data(), aa_seq_len);
    PrintNucCodeSequence(nuc_seq_str.get(), nuc_sequence.data(), nuc_seq_len);
    
    int n_evaluations =
	n_gens_without_accept * (int)json_log["reevaluation size"]
	+ (int)json_log["simulations per gen"];
    if (generation == 0)
    {
	n_evaluations = 0;
    }

    json entry;
    entry["generation"] = generation;
    entry["old evaluations"] = n_evaluations;
    entry["mutation type"] = mutation_type;
    entry["aa sequence"] = aa_seq_str.get();
    entry["nuc sequence"] = nuc_seq_str.get();
    entry["old fitness"] = old_fitness;
    entry["new fitness"] = new_fitness;
    entry["old protein output"] = old_protein_output;
    entry["new protein output"] = new_protein_output;
    entry["old pnat"] = old_pnat;
    entry["new pnat"] = new_pnat;
    entry["native energy"] = native_energy;
    entry["accepted"] = accept;
    json_log["trajectory"].push_back(entry);
    return;
}


// Write json file (does overwrite)
void write_log(
    const json& json_log,
    const std::string& json_log_path)
{
    if (!g_world_rank)
    {
	std::fstream out(json_log_path, std::ios::out | std::ios::trunc);
	if (!out.is_open())
	{
	    std::cerr << "Output file could not be opened";
	    exit(1);
	}

	out << std::setw(2) << json_log << std::endl;
	out.close();
    }
    return;
}
