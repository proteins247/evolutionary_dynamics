/* 
 * folding_evolution v0.0.4
 *
 * v0.0.4 implements checkpointing
 * 
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

#include <getopt.h>
#include <unistd.h>
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
    "  -s, --speed-params=FILE   REQUIRED. Path to file with translation\n"
    "                            times of the codons.\n"
    "  -c, --save-conformations  Save latFoldVec simulation conformations.\n"
    "  -d, --debug               Allow error messages from all processors.\n"
    "  -f, --from-checkpoint     Resume simulation from checkpoint file.\n"
    "  -k, --degradation-param=K Value setting timescale for unfolded protein\n"
    "                            degradation (penalizes slow folding).\n"
    "  -l, --latpack-path=PATH   Latpack binaries here. Default=\n"
    "                            /n/home00/vzhao/pkg/latPack/1.9.1-8/\n"
    "  -m, --mutation-mode=M     One of 0, 1, or 2, indicating restriction to\n"
    "                            synonymous mutation, nonsynonymous mutation,\n"
    "                            or allow both kinds of mutations (default).\n"
    "  -o, --out-path=FILE       Output will be put in this directory.\n"
    "                            Default=./out\n"
    "  -p, --population-size=N   Size of population for evolutionary\n"
    "                            dynamics. Default=500\n"
    "      --random-codons       If SEQUENCE is of amino acids, corresponding\n"
    "                            RNA sequence codons will be randomly chosen.\n"
    "                            (Default behavior is to use fastest codons.)\n"
    "  -r, --seed=N              RNG seed. Default=1\n"
    "  -t, --temperature=T       Temperature of latFoldVec simulations.\n"
    "                            Default=0.3\n"
    "      --help   Display this help and exit.\n"
    "\n"
    "Format of the output files:\n"
    "\nThe primary output is a simulation log that records the generation\n"
    "number, protein sequence, nucleic acid sequence, fitness, and whether\n"
    "the mutation was accepted. Trajectory files are also saved.\n"
    "\n"
    ;

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;
static const int LATPACK_ERROR = 4;

// default values
static const int GENS_MAX = 99999;
static const std::string DEFAULT_LATPACK_PATH =
    "/n/home00/vzhao/pkg/latPack/1.9.1-9/";
static const std::string DEFAULT_OUTPATH = "./out";
static const uint64_t DEFAULT_SEED = 1;
static const int DEFAULT_LATFOLD_OUTFREQ = 10000;
static const int DEFAULT_CHECKPOINT_FREQ = 20;
static const int DEFAULT_JSON_OUTFREQ = 5;
static const int DEFAULT_POPULATION_SIZE = 500;
static const int DEFAULT_MUTATION_MODE = 2;    // MutateAll default value
static const double DEFAULT_TEMPERATURE = 0.3;
static const double DEFAULT_REEVALUATION_RATIO = 0.5; // Not commandline option.
static const double DEFAULT_FITNESS_CONSTANT = 0.75;
static const double DEFAULT_DEGRADATION_PARAM = 1000000;

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


// Given a codon sequence and a dictionary linking codons and
// translation times, make a string of space-separated translation times.
//
// Codons codon_sequence[5] through
// codon_sequence[codon_sequence.size()] will be used, and the last
// time used will be final_time.
// 
// @param times Codon to translation time dictionary
// @param codon_sequence The sequence of codons
// @param final_time The last number to put in the string (representing
//        time to release from ribosome)
// @param total_translation_steps This function places the number of
//        steps translation will take into this variable.
// 
// @return A string of space-separated translation times corresponding
//   to the given codon sequence to be used for the
//   -elongationSchedule argument in latFoldVec.
std::string make_translation_schedule(
    const std::vector<unsigned int> & times,
    const std::vector<Codon> & codon_sequence,
    const int final_time,
    int* total_translation_steps);


// Build a vector of strings that are arguments to invoke
// latFoldVec. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfoldvec_command(
    const std::string& latpack_path,
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    const std::string& translation_schedule,
    const int& full_length_time,
    const double& temperature,
    const int& output_frequency,
    bool save_conformations);


// Use vfork and exec to run a latFoldVec co-translational folding
// simulation.
//
// @param latfoldvec_command The vector of strings built
//        by compose_latfoldvec_command.
// @param outfile latFoldVec program output will go to this file.
void run_latfoldvec(
    std::vector<std::string> latfoldvec_command,
    const std::string & outfile);


// Build a vector of strings that are arguments to invoke
// latMapTraj. Note arguments -traj is not added here.
std::vector<std::string> compose_latmaptraj_command(
    const std::string& latpack_path,
    const std::string& folded_conformation);


// FIXME (this function is not currently used)
// Use fork and exec to run n_simulations of latMapTraj to analyze the
// latFoldVec simulations that were run.
void run_latmaptraj(
    std::vector<std::string> & latmaptraj_command,
    const std::string& h5file_base,
    int n_simulations=1);


// Analyze latFoldVec simulations to count number of successful
// simulations. This is a parallel function.
//
// @param filename The partial name of the file hdf5 output was
//        written to.
// @param total_translation_steps Total number of MC steps needed
// @param degradation_param Timescale that unfolded protein degradation
//        acts on.
// @param native_energy Value is updated with native energy.
double evaluate_folded_fraction(
    const std::string& filename,
    int total_translation_steps,
    double degradation_param,
    double* native_energy);


// Calculate fitness from folded fraction.
//
// @param folded_fraction The fraction of successful folding
//        simulations
// @param f_0 Fitness function parameter
double calculate_fitness(
    double folded_fraction,
    double f_0=DEFAULT_FITNESS_CONSTANT);


// Update folded fraction to incorporate new folding simulations.
//
// @param old_folded_fraction The previously estimated folded fraction.
// @param new_folded_fraction The newly estimated folded fraction.
// @param n_reevaluators The number of simulations run to estimate
//        the new fitness value.
// @param n_total The total number of CPUs running simulations.
// @param n_gens_without_mutation Number of generations since last
//        accepted mutation.
double reaverage_folded_fraction(
    double old_folded_fraction,
    double new_folded_fraction,
    int n_reevaluators,
    int n_total,
    int n_gens_without_mutation);


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
    int n_gens_without_mutation,
    int mutation_type,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double old_folded_fraction,
    double new_folded_fraction,
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
    bool proc_does_reevaluation = (g_world_rank < reevaluation_size) ?
	true : false;
    MPI_Comm_split(
	MPI_COMM_WORLD, proc_does_reevaluation, g_world_rank, &g_subcomm);
    MPI_Comm_rank(g_subcomm, &g_subcomm_rank);
    MPI_Comm_size(g_subcomm, &g_subcomm_size);
    assert(g_subcomm_size > 0);

    // Variables to be determined by commandline options: 

    // Required:
    std::string sequence;
    int n_gens;			// How many generations to run

    // Other params:
    bool debug_mode = false;
    bool save_conformations = false;
    bool resume_from_checkpoint = false;
    int random_codons = 0;
    int population_size = DEFAULT_POPULATION_SIZE;
    std::string folded_conformation;
    std::string speedparams_path;
    std::string latpack_path = DEFAULT_LATPACK_PATH;
    std::string out_path = DEFAULT_OUTPATH;
    std::string checkpoint_path;
    std::string lat_sim_out_path;
    std::string json_log_path;	// we're going to log using json
    json checkpoint;
    json json_log;
    std::ostream * outstream = &std::cout;
    double degradation_param = DEFAULT_DEGRADATION_PARAM;
    uint64_t rng_seed = DEFAULT_SEED;
    enum MutationMode
    {
	SynonymousOnly,
	NonsynonymousOnly,
	MutateAll
    } mutation_mode = static_cast<MutationMode>(DEFAULT_MUTATION_MODE);
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
	{"from-checkpoint", required_argument, NULL, 'f'},
	{"degradation-param", required_argument, NULL, 'k'},
	{"latpack-path", required_argument, NULL, 'l'},
	{"mutation-mode", required_argument, NULL, 'm'},
	{"native-fold", required_argument, NULL, 'n'},
	{"random-codons", no_argument, &random_codons, 1},
	{"out-path", required_argument, NULL, 'o'},
	{"population-size", required_argument, NULL, 'p'},
	{"seed", required_argument, NULL, 'r'},
	{"speed-params", required_argument, NULL, 's'},
	{"temperature", required_argument, NULL, 't'},
	{"save-conformations", no_argument, NULL, 'c'},
	{"debug", no_argument, NULL, 'd'},
	{"help", no_argument, NULL, 'h'},
	{NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "f:k:l:m:n:o:p:r:s:t:hdc",
			    long_options, &option_index))
	   != -1)
    {
	switch (c)
	{
	case 'h':
	    print_help();
	    exit(0);
	case 'f':
	    // We are resuming from checkpoint. Ignore all other arguments
	    checkpoint_path = optarg;
	    resume_from_checkpoint = true;
	    goto checkpoint_breakout;
	case 'c':
	    save_conformations = true;
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
	case 'l':
	    latpack_path = optarg;
	    break;
	case 'm':
	    try
	    {
		switch (std::stoi(optarg))
		{
		case 0:
		    mutation_mode = SynonymousOnly;
		    break;
		case 1:
		    mutation_mode = NonsynonymousOnly;
		    break;
		case 2:
		    mutation_mode = MutateAll;
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
	case 's':
	    speedparams_path = optarg;
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
	continue;
    checkpoint_breakout:
	break;
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
	if (speedparams_path.empty())
	{
	    std::cerr << "Argument --speed-params (-s) is mandatory"
		      << std::endl;
	    exit(PARSE_ERROR);
	}

	lat_sim_out_path = out_path + "/latfoldvec_simulations";
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

	// Wait for slow network file system
	MPI_Barrier(MPI_COMM_WORLD);
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
	json_log.at("generations").get_to(n_gens);
	json_log.at("population size").get_to(population_size);
	json_log.at("temperature").get_to(temperature);
	json_log.at("simulations per gen").get_to(simulations_per_gen);
	json_log.at("reevaluation size").get_to(checkpoint_reevaluation_size);
	json_log.at("degradation timescale").get_to(degradation_param);
	json_log.at("latpack path").get_to(latpack_path);
	json_log.at("translation params").get_to(speedparams_path);
	json_log.at("mutation mode").get_to(mutation_mode);

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
	int checkpoint_generation;
	checkpoint.at("generation").get_to(checkpoint_generation);
	auto it = json_log.at("trajectory").begin();
	while (it != json_log.at("trajectory").end() &&
	       it->at("generation") <= checkpoint_generation)
	{
	    ++it;
	}
	json_log.at("trajectory").erase(it, json_log.at("trajectory").end());

	// for (auto it = json_log.at("trajectory").begin();
	//      it != json_log.at("trajectory").end();
	//      ++it)
	// {
	//     if (it->at("generation") > checkpoint_generation)
	//     {
	// 	json_log.at("trajectory").erase(it);
	//     }
	// }
    }
    
    // End option parsing
    // --------------------------------------------------
    // Begin simulation setup

    // Process the user-provided sequence
    std::vector<int> nuc_sequence;
    std::vector<int> prev_nuc_sequence;
    std::vector<AminoAcid> aa_sequence;
    std::vector<Codon> codon_sequence;

    // Infer lengths from folded conformation length
    unsigned int protein_length = folded_conformation.length() + 1;
    unsigned int nuc_length = protein_length * 3;

    // Resize vectors
    nuc_sequence.resize(nuc_length);
    prev_nuc_sequence.resize(nuc_length);
    aa_sequence.resize(protein_length);
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
	prev_nuc_sequence.clear();
	checkpoint.at("nuc sequence").get_to(nuc_sequence);
	checkpoint.at("previous nuc sequence").get_to(prev_nuc_sequence);
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());
    }

    // Obtain a codon sequence from nuc_sequence
    NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());

    // Open translation speed data file and read
    // Make vector of 821, since highest codon is 0x333 = 820
    std::vector<unsigned int> translation_times(821, 0);
    std::vector<int> stop_codon_times;
    ReadTranslationTimes(speedparams_path.c_str(), translation_times.data());

    // Here, we find the three translation times that correspond to
    // the stop codons.
    for (auto& it : STOP_CODONS)
    {
	stop_codon_times.push_back(translation_times.at(it));
    }

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
	    json_log["latpack path"] = latpack_path;
	    json_log["translation params"] = speedparams_path;
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
	    << "# rng seed : " << rng_seed << std::endl
	    << "# translation params : " << speedparams_path << std::endl
	    ;
	print_header(outstream, aa_sequence, nuc_sequence); 
    }
    
    // Simulation-related variables and parameters
    // First determine the time allowed for posttranslational folding.
    // The probability of posttranslation degradation is 1 - exp(-t/tau)
    // where tau = degradation_param and t = time since translation
    // (ribosome release).
    // Set it to be log(8) * tau. (87.5% probability of degradation).
    int posttranslational_folding_time = log(8) * degradation_param;
    int gen = 0;
    int old_total_translation_steps;
    int total_translation_steps;
    int last_accepted_gen = 0;
    int n_gens_without_mutation = 0;
    int mutation_type = -1;
    double old_fitness = 0.001;
    double fitness;
    double old_folded_fraction = 0;
    double folded_fraction;
    double native_energy;
    double selection;
    double fixation;
    std::string translation_schedule;
    std::vector<std::string> latfoldvec_command;
    std::vector<std::string> prev_latfoldvec_command;
    std::unique_ptr<char []> aa_sequence_str(new char[protein_length+1]);

    if (resume_from_checkpoint)
    {
	checkpoint.at("generation").get_to(gen);
	++gen;
	checkpoint.at("old total translation steps").get_to(
	    old_total_translation_steps);
	checkpoint.at("last accepted gen").get_to(last_accepted_gen);
	checkpoint.at("n gens without mutation").get_to(
	    n_gens_without_mutation);
	checkpoint.at("mutation type").get_to(mutation_type);
	checkpoint.at("old fitness").get_to(old_fitness);
	checkpoint.at("old folded fraction").get_to(old_folded_fraction);
	checkpoint.at("prev latfoldvec command").get_to(
	    prev_latfoldvec_command);
    }

    // Begin running simulation loop
    for (; gen < n_gens; ++gen)
    {
	// Compose the simulation parameters.
	translation_schedule = make_translation_schedule(
	    translation_times, codon_sequence, stop_codon_times[0],
	    &total_translation_steps);
	PrintAASequence(
	    aa_sequence_str.get(), aa_sequence.data(), protein_length);
	latfoldvec_command = compose_latfoldvec_command(
	    latpack_path,
	    aa_sequence_str.get(),
	    folded_conformation,
	    translation_schedule,
	    posttranslational_folding_time,
	    temperature,
	    latfold_output_frequency,
	    save_conformations);

	if (proc_does_reevaluation && gen > 0)
	{
	    std::ostringstream output_dir;
	    std::ostringstream hdf5_output_file;
	    output_dir << lat_sim_out_path << "/gen" << std::setw(5)
		       << std::setfill('0') << last_accepted_gen;
	    hdf5_output_file << output_dir.str() << "/sim" << std::setfill('0')
			     << std::setw(5) << n_gens_without_mutation + 1
			     << "_" << std::setw(5) << g_subcomm_rank << ".h5";

	    run_latfoldvec(prev_latfoldvec_command, hdf5_output_file.str());
	    folded_fraction = evaluate_folded_fraction(
		hdf5_output_file.str(), old_total_translation_steps,
		degradation_param, &native_energy);

	    // Now update old fitness
	    old_folded_fraction = reaverage_folded_fraction(
		old_folded_fraction, folded_fraction, g_subcomm_size,
		g_world_size, n_gens_without_mutation);
	    old_fitness = calculate_fitness(old_folded_fraction);
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

	    run_latfoldvec(latfoldvec_command, hdf5_output_file.str());
	    folded_fraction = evaluate_folded_fraction(
		hdf5_output_file.str(), total_translation_steps,
		degradation_param, &native_energy);
	    fitness = calculate_fitness(folded_fraction);
	}
	MPI_Bcast(&native_energy, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&fitness, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&folded_fraction, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&old_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&old_folded_fraction, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
	    save_state(json_log, gen, n_gens_without_mutation, mutation_type,
		       aa_sequence, nuc_sequence, old_fitness, fitness,
		       old_folded_fraction, folded_fraction, native_energy,
		       (bool)accept);
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
	    // Update so that nuc_sequence and fitness are fixed
	    prev_nuc_sequence = nuc_sequence;
	    prev_latfoldvec_command = latfoldvec_command;
	    old_total_translation_steps = total_translation_steps;
	    old_folded_fraction = folded_fraction;
	    old_fitness = fitness;
	    n_gens_without_mutation = 0;
	    last_accepted_gen = gen;
	}
	else
	{
	    // Revert
	    nuc_sequence = prev_nuc_sequence;
	    n_gens_without_mutation++;
	}

	// Make new mutation. We only do this on root node because of
	// independent rng streams for different nodes.
	if (!g_world_rank)
	{
	    do
	    {
		std::vector<int> temp_sequence = nuc_sequence;
		mutation_type = PointMutateNucSequence(
		    nuc_sequence.data(), nuc_length);
		if (mutation_type == mutation_mode)
		    break;
		if (mutation_mode == MutateAll && mutation_type >= 0)
		    break;
		nuc_sequence = temp_sequence;
	    } while (true);
	}

	MPI_Bcast(nuc_sequence.data(), nuc_sequence.size(), MPI_INT,
		  0, MPI_COMM_WORLD);
	
	// and update sequences
	NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());

	// Need to save checkpoint
	if (gen % checkpoint_frequency == 0)
	{
	    std::ostringstream checkpoint_path;
	    checkpoint_path << out_path << "/checkpoint"<< std::setw(5)
			   << std::setfill('0') << std::to_string(gen)
			   << ".json";

	    checkpoint["generation"] = gen;
	    checkpoint["nuc sequence"] = nuc_sequence;
	    checkpoint["previous nuc sequence"] = prev_nuc_sequence;
	    checkpoint["old total translation steps"] =
		old_total_translation_steps;
	    checkpoint["last accepted gen"] = last_accepted_gen;
	    checkpoint["n gens without mutation"] =
		n_gens_without_mutation;
	    checkpoint["mutation type"] = mutation_type;
	    checkpoint["old fitness"] = old_fitness;
	    checkpoint["old folded fraction"]
		= old_folded_fraction;
	    checkpoint["prev latfoldvec command"] = prev_latfoldvec_command;
	    write_checkpoint(checkpoint_path.str(), checkpoint);
	}
    }

    // finalize json log
    write_log(json_log, json_log_path);
    
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


// Given a codon sequence and a dictionary linking codons and
// translation times, make a string of space-separated translation times.
std::string make_translation_schedule(
    const std::vector<unsigned int> & times,
    const std::vector<Codon> & codon_sequence,
    const int final_time,
    int* total_translation_steps)
{
    std::ostringstream os;
    *total_translation_steps = 0;
    int protein_length = 5;

    // Note that we skip the first five codons (start at the
    //   6th--index 5), since those residues are present at the start
    //   of the MC simulation.
    for (auto codon = codon_sequence.begin() + 5;
	 codon != codon_sequence.end(); ++codon)
    {
	os << times.at(*codon) << " ";
	*total_translation_steps += times.at(*codon) * protein_length;
	++protein_length;
    }
    os << final_time;
    *total_translation_steps += final_time * protein_length;
    return os.str();
}


// Build a vector of strings that are arguments to invoke
// latFoldVec. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfoldvec_command(
    const std::string& latpack_path,
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    const std::string& translation_schedule,
    const int& full_length_time,
    const double& temperature,
    const int& output_frequency,
    bool save_conformations)
{
    std::vector<std::string> command;

    command.push_back(latpack_path + "/bin/latFoldVec");
    command.push_back("-energyFile=" + latpack_path + "/share/latpack/MJ.txt");
    command.push_back("-seq=" + aa_sequence);
    command.push_back("-final=" + folded_conformation);
    command.push_back("-elongationSchedule=" + translation_schedule);
    command.push_back("-kT=" + std::to_string(temperature));
    command.push_back("-maxStepsIncrease");
    command.push_back("-ribosome");
    command.push_back("-ribosomeRelease");
    command.push_back("-fullLengthSteps=" + std::to_string(full_length_time));
    command.push_back("-outFreq=" + std::to_string(output_frequency));
    if (save_conformations)
    {
	command.push_back("-out=S");
    }
    else
    {
	command.push_back("-out=E");
    }
    
    // Previously we did the what is commented out, but we want to
    // keep the command as a vector.
    // std::ostringstream os;
    // std::copy(command.begin(), command.end(),
    // std::ostream_iterator<std::string>(os, " ")); return os.str();
    return command;
}


// Use fork and exec to run a latFoldVec co-translational folding
// simulation.
void run_latfoldvec(
    std::vector<std::string> latfoldvec_command,
    const std::string & outfile)
{
    // Add additional parameters
    latfoldvec_command.push_back(
	"-seed=" + std::to_string(threefryrand_int() % INT32_MAX));
    latfoldvec_command.push_back("-outFile=" + outfile);
    // latfoldvec_command.push_back("-title=?");

    std::ostringstream command;
    for (auto& arg : latfoldvec_command)
    {
	command << arg << " ";
    }

    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
	latfoldvec_command);

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

	execv(cstring_command_vec[0], cstring_command_vec.data());
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
// latFoldVec simulations that were run.
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



// Analyze latFoldVec simulations to determine a fitness value for
// protein under evaluation.
double evaluate_folded_fraction(
    const std::string& filename,
    int total_translation_steps,
    double degradation_param,
    double* native_energy)
{
    int folded = 0;
    int total_n_folded = 0;

    // We need to read data from HDF5 file
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_id = H5Gopen2(file_id, "traj1", H5P_DEFAULT);
    hid_t sequence_attr = H5Aopen(file_id, "Sequence", H5P_DEFAULT);
    hid_t sequence_attr_t = H5Aget_type(sequence_attr);
    size_t sequence_attr_size = H5Tget_size(sequence_attr_t);
    size_t protein_length = sequence_attr_size - 1;
    hid_t found_final_attr = H5Aopen(group_id, "Found final struct",
				     H5P_DEFAULT);
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 4);
    hid_t last_step_attr = H5Aopen(group_id, "Last step", H5P_DEFAULT);
    hid_t last_energy_attr = H5Aopen(group_id, "Last E", H5P_DEFAULT);

    char found_final[4];
    int last_step;
    double last_energy;

    H5Aread(found_final_attr, strtype, found_final);
    H5Aread(last_energy_attr, H5T_NATIVE_DOUBLE, &last_energy);
    H5Aread(last_step_attr, H5T_NATIVE_UINT, &last_step);

    if (std::strncmp(found_final, "yes", 3) == 0)
    {
	int posttranslation_steps = last_step - total_translation_steps;
	if (posttranslation_steps < 0)
	{
	    folded = 1;
	}
	else
	{
	    double posttranslation_time =
		posttranslation_steps / protein_length;
	    double degradation_probability = 1 - exp(
		-posttranslation_time / degradation_param);
	    if (threefryrand() > degradation_probability)
	    {
		folded = 1;
	    }
	}
    }
    else
    {
	last_energy = 0;
    }
	
    H5Tclose(strtype);
    H5Tclose(sequence_attr_t);
    H5Aclose(sequence_attr);
    H5Aclose(found_final_attr);
    H5Aclose(last_step_attr);
    H5Aclose(last_energy_attr);    
    H5Gclose(group_id);
    H5Fclose(file_id);

    if (g_subcomm_rank)
    {
	// Send whether trajectory folded
	MPI_Send(&folded, 1, MPI_INT, 0, 0, g_subcomm);
    }
    else  // root node
    {
	total_n_folded += folded;
	for (int i=1; i<g_subcomm_size; ++i)
	{
	    MPI_Recv(&folded, 1, MPI_INT, MPI_ANY_SOURCE, 0,
		     g_subcomm, &g_status);
	    total_n_folded += folded;
	}
    }
    // Let everyone know total folded
    MPI_Bcast(&total_n_folded, 1, MPI_INT, 0, g_subcomm);
    double folded_fraction = total_n_folded / (double)g_subcomm_size;

    // Also determine native energy
    MPI_Allreduce(&last_energy, native_energy, 1, MPI_DOUBLE, MPI_MIN,
		  g_subcomm);
    
    return folded_fraction;
}


// Calculate fitness from folded fraction.
double calculate_fitness(
    double folded_fraction,
    double f_0)
{
    return 0.001 + folded_fraction / (folded_fraction + f_0);    
}


// Update folded fraction to incorporate new folding simulations.
double reaverage_folded_fraction(
    double old_folded_fraction,
    double new_folded_fraction,
    int n_reevaluators,
    int n_total,
    int n_gens_without_mutation)
{
    // We calculate a weighted average, so here we get the weights.
    double old_factor = (n_total + n_reevaluators * (n_gens_without_mutation - 1))
	/ (double)(n_total + n_gens_without_mutation * n_reevaluators);
    double new_factor = n_reevaluators
	/ (double)(n_total + n_gens_without_mutation * n_reevaluators);
    return old_folded_fraction * old_factor + new_folded_fraction * new_factor;
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

    *outstream << std::fixed << std::setprecision(2);
    
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
	       << std::setw(8) << old_fitness << " "
	       << std::setw(8) << new_fitness << " "
	       << std::setw(7) << native_energy << " "	
	       << std::setw(6) << ((accept) ? "yes" : "no")
	       << std::endl;

    return;    
}


// Save results to json log.
void save_state(
    json& json_log,
    int generation,
    int n_gens_without_mutation,
    int mutation_type,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    double old_folded_fraction,
    double new_folded_fraction,
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
	n_gens_without_mutation * (int)json_log["reevaluation size"]
	+ ((int)json_log["simulations per gen"]
	   - (int)json_log["reevaluation size"]);

    json entry;
    entry["generation"] = generation;
    entry["old evaluations"] = n_evaluations;
    entry["mutation type"] = mutation_type;
    entry["aa sequence"] = aa_seq_str.get();
    entry["nuc sequence"] = nuc_seq_str.get();
    entry["old fitness"] = old_fitness;
    entry["new fitness"] = new_fitness;
    entry["old folded fraction"] = old_folded_fraction;
    entry["new folded fraction"] = new_folded_fraction;
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
