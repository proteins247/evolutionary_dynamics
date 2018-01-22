/* 
 * folding_evolution
 *
 * Monoclonal simulation
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
#include <cmath>
#include <map>

#include <getopt.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include <hdf5.h>
#include <mpi.h>

#define CPLUSPLUS __cplusplus
#undef __cplusplus
// When compiling Random123 header files for rng.h, we don't need the
//   C++ features (and they're not compatible with extern "C"), so we
//   temporarily undefine __cplusplus, saving its value in CPLUSPLUS.
//   Alternatively, we could build rng.o from rng.c using g++, and then
//   we wouldn't have to wrap rng.h in extern "C".

extern "C"
{
#include "../gencode.h"
#include "../rng.h"
}

#define __cplusplus CPLUSPLUS
#undef CPLUSPLUS


static const std::string helptext =
    "folding_evolution\n"
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
    "generation, the previous sequence is also reevaluated.\n"
    "\n"
    "  -n, --native-fold=CONF    Required. The native conformation of the\n"
    "                            as a latPack-style move sequence.\n"
    "  -s, --speed-params=FILE   Required. Path to file with translation\n"
    "                            times of the codons.\n"
    "  -d, --h5file-path=PATH    LatFoldVec trajectory files written to this\n"
    "                            directory. Directory must exist. Default=./\n"
    "  -l, --latpack-path=PATH   Latpack binaries here. Default=\n"
    "                            /n/home00/vzhao/pkg/latPack/1.9.1-7/\n"
    "  -m, --mutation-mode=M     One of 0, 1, or 2, indicating restriction to\n"
    "                            synonymous mutation, nonsynonymous mutation,\n"
    "                            or allow both kinds of mutations (default).\n"
    "  -o, --outfile=FILE        Output will be sent to this file instead of\n"
    "                            STDOUT.\n"
    "  -p, --population-size=N   Size of population for evolutionary\n"
    "                            dynamics. Default=500\n"
    "      --random-codons       If SEQUENCE is of amino acids, corresponding\n"
    "                            RNA sequence codons will be randomly chosen.\n"
    "                            (Default behavior is using fastest codons.)\n"
    "  -r, --seed=N              RNG seed. Default=1\n"
    "  -t, --temperature=T       Temperature of latFoldVec simulations.\n"
    "                            Default=0.3\n"
    "      --help   Display this help and exit.\n"
    "\n"
    "Format of the output files:\n"
    "\nThe primary output is a simulation log that records the generation\n"
    "number, protein sequence, nucleic acid sequence, fitness, and whether\n"
    "the mutation was accepted.\n"
    "\n"
    "In addition, it's recommended that --h5file-path be set so that\n"
    "trajectory files go into a single directory. (The default behavior will\n"
    "put files into the user's current directory.)\n"
    "\n"
    ;

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;
static const int LATPACK_ERROR = 4;

// default values
static const std::string DEFAULT_LATPACK_PATH = "/n/home00/vzhao/pkg/latPack/1.9.1-7/";
static const std::string DEFAULT_H5FILE_PATH = "./";
static const uint64_t DEFAULT_SEED = 1;
static const int DEFAULT_LATFOLD_OUTFREQ = 1000;
static const int DEFAULT_POPULATION_SIZE = 500;
static const int DEFAULT_MUTATION_MODE = 2;    // MutateAll default value
static const double DEFAULT_TEMPERATURE = 0.3;
static const double DEFAULT_REEVALUATION_RATIO = 0.5;

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
// @param final_time The last number to put in the string
// 
// @return A string of space-separated translation times corresponding
//   to the given codon sequence to be used for the
//   -elongationSchedule argument in latFoldVec.
std::string make_translation_schedule(
    const std::vector<unsigned int> & times,
    const std::vector<Codon> & codon_sequence,
    const int & final_time=0);


// Build a vector of strings that are arguments to invoke
// latFoldVec. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfoldvec_command(
    const std::string& latpack_path,
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    const std::string& translation_schedule,
    const double& temperature,
    const int& output_frequency);


// Use fork and exec to run a latFoldVec co-translational folding
// simulation.
//
// @param latfoldvec_command The vector of strings built
//        by compose_latfoldvec_command.
// @param h5file_base latFoldVec program output will go to
//        a file whose pathname begins with this.
void run_latfoldvec(
    std::vector<std::string> & latfoldvec_command,
    const std::string & h5file_base);


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
// @param h5file_base The partial name of the file hdf5 output was
//        written to.
// @param n_0 Fitness function parameter.
double evaluate_folded_fraction(
    const std::string& filename);


// Calculate fitness from folded fraction.
//
// @param folded_fraction The fraction of successful folding
//        simulations
// @param f_0 Fitness function parameter
double calculate_fitness(
    double folded_fraction,
    double f_0=0.5);


// Update fitness to incorporate new fitness evaluations.
//
// @param old_fitness The previously estimated fitness value.
// @param new_fitness The newly estimated fitness value.
// @param n_reevaluators The number of simulations run to estimate
//        the nw fitness value.
// @param n_total The total number of CPUs running simulations.
// @param n_gens_without_mutation Number of generations since last
//        accepted mutation.
double reaverage_folded_fraction(
    double old_folded_fraction,
    double new_folded_fraction,
    int n_reevaluators,
    int n_total,
    int n_gens_without_mutation);


// Print the header for state output
// we would want
void cout_header(
    std::ostream * outstream,
    const std::vector<AminoAcid> & aa_sequence,
    const std::vector<int> & nuc_sequence);


// Print out what happened each generation.
// todo is to handle sequence reevaluation
void cout_state(
    std::ostream * outstream,
    int generation,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
    bool accept);


int main(int argc, char** argv)
{
    // MPI init stuff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &g_world_size);
    assert (g_world_size > 1);    // We need at least two CPUs

    // Split cpus into two groups
    // world_rank [0, reevaluation_size) does reevaluation
    int reevaluation_size = (int)g_world_size * DEFAULT_REEVALUATION_RATIO;
    bool proc_does_reevaluation = (g_world_rank < reevaluation_size) ?
	true : false;
    MPI_Comm_split(
	MPI_COMM_WORLD, proc_does_reevaluation, g_world_rank, &g_subcomm);
    MPI_Comm_rank(g_subcomm, &g_subcomm_rank);
    MPI_Comm_size(g_subcomm, &g_subcomm_size);

    // Variables to be determined by commandline options: 
    std::string sequence;
    int n_gens;			// How many generations to run
    std::string folded_conformation;
    std::string speedparams_path;
    std::string outfile_path;
    std::string latpack_path = DEFAULT_LATPACK_PATH;
    std::string h5file_path = DEFAULT_H5FILE_PATH;
    std::ostream * outstream = &std::cout;
    std::ofstream outfile;
    uint64_t rng_seed = DEFAULT_SEED;
    enum MutationMode
    {
	SynonymousOnly,
	NonsynonymousOnly,
	MutateAll
    } mutation_mode = static_cast<MutationMode>(DEFAULT_MUTATION_MODE);
    double temperature = DEFAULT_TEMPERATURE;

    // These particular params currently cannot be changed by
    // commandline arguments
    int latfold_output_frequency = DEFAULT_LATFOLD_OUTFREQ;
    int random_codons = 0;
    int population_size = DEFAULT_POPULATION_SIZE;
    // End variables to be determined by commandline options: 

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
	{"h5file-path", required_argument, NULL, 'd'},
	{"latpack-path", required_argument, NULL, 'l'},
	{"mutation-mode", required_argument, NULL, 'm'},
	{"native-fold", required_argument, NULL, 'n'},
	{"random-codons", no_argument, &random_codons, 1},
	{"outfile", required_argument, NULL, 'o'},
	{"population-size", required_argument, NULL, 'p'},
	{"seed", required_argument, NULL, 'r'},
	{"speed-params", required_argument, NULL, 's'},
	{"temperature", required_argument, NULL, 't'},
	{"help", no_argument, NULL, 'h'},
	{NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "d:l:m:n:o:p:r:s:t:h",
			    long_options, &option_index))
	   != -1)
    {
	switch (c)
	{
	case 'h':
	    print_help();
	    exit(0);
	case 0:
	    break;
	case 'd':
	    h5file_path = optarg;
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
		std::cerr << "--mutation-mode argument not one of 0, 1, 2: "
			  << optarg << std::endl;
		exit(PARSE_ERROR);
	    }
	    break;
	case 'n':
	    folded_conformation = optarg;
	    break;
	case 'o':
	    outfile_path = optarg;
	    break;
	case 'p':
	    try
	    {
		population_size = std::stoi(optarg);
	    }
	    catch (...)
	    {
		std::cerr << "Failed to convert --population-size: " << optarg
			  << std::endl;
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
		std::cerr << "Failed to convert seed: " << optarg << std::endl;
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
		std::cerr << "Failed to convert temp: " << optarg << std::endl;
		exit(PARSE_ERROR);
	    }
	    break;
	case '?':
	    // getopt_long will print an error message
	    print_help();
	    exit(PARSE_ERROR);
	default:
	    // we shouldn't end up here i think
	    std::cerr << "getopt switch default; returned " << c << std::endl;
	    exit(PARSE_ERROR);
	}
    } // End while getopt_long

    // Continued parsing: positional arguments
    if (optind + 2 == argc)
    {
	sequence = argv[optind++];
	try
	{
	    n_gens = std::stoi(argv[optind]);
	}
	catch (const std::invalid_argument& e)
	{
	    std::cerr << "Provided N_GENS argument " << argv[optind]
		      << " could not be converted to uint" << std::endl;
	    exit(PARSE_ERROR);
	}
	catch (const std::out_of_range& e)
	{
	    // probably unlikely
	    std::cerr << "Provided N_GENS argument " << argv[optind]
		      << " is out of range" << std::endl;
	    exit(PARSE_ERROR);
	}
    }
    else
    {
	std::cerr << "Expected two positional arguments: SEQUENCE NROUNDS"
		  << std::endl;
	print_help();
	exit(PARSE_ERROR);
    }
    // Complain about missing parameters here
    if (folded_conformation.empty())
    {
	std::cerr << "Argument --native-fold (-n) is mandatory" << std::endl;
	exit(PARSE_ERROR);
    }
    if (speedparams_path.empty())
    {
	std::cerr << "Argument --speed-params (-s) is mandatory" << std::endl;
	exit(PARSE_ERROR);
    }
    
    // End option parsing
    // --------------------------------------------------
    // Begin simulation setup

    // Initialize RNG, giving different processes different seeds
    set_threefry_array(rng_seed, g_world_rank, proc_does_reevaluation, 0);

    // Setup output
    if (!g_world_rank && !outfile_path.empty())
    {
	outfile.open(outfile_path);
	if (!outfile)
	{
	    std::cerr << "Problem with output file: "
		      << outfile_path << std::endl;
	    exit(IO_ERROR);
	}
	outstream = &outfile;
    }
    
    // Process the user-provided sequence
    std::vector<int> nuc_sequence;
    std::vector<int> prev_nuc_sequence;
    std::vector<AminoAcid> aa_sequence;
    std::vector<Codon> codon_sequence;

    // Infer lengths from folded conformation length
    unsigned int protein_length = folded_conformation.length() + 1;
    unsigned int nuc_length = protein_length * 3;

    // After inferring length, resize vectors
    nuc_sequence.resize(nuc_length);
    prev_nuc_sequence.resize(nuc_length);
    aa_sequence.resize(protein_length);
    codon_sequence.resize(protein_length);

    // We use resize because we directly access the underlying arrays of
    // our three vectors in this next section
    if (protein_length == sequence.length()) 
    {
	// presume it's a protein sequence
	if (LetterToAASeq(
		sequence.c_str(), aa_sequence.data(), protein_length))
	{
	    std::cerr << "Letter to aa conversion problem" << std::endl;
	    std::cerr << "Sequence: " << sequence << std::endl;
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
	    std::cerr << "letter to nuc conversion problem" << std::endl;
	    std::cerr << "Sequence: " << sequence << std::endl;
	    exit(DATA_ERROR);
	}
	if (NucSeqToAASeq(
		nuc_sequence.data(), nuc_length, aa_sequence.data()))
	{
	    // nuc sequence contains a stop codon
	    std::cerr << "Provided RNA sequence contains a stop codon"
		      << std::endl;
	    std::cerr << "Sequence: " << sequence << std::endl;
	    exit(DATA_ERROR);
	}
    }
    else
    {
	// houston we have a problem
	std::cerr << "Protein length of " << protein_length
		  << ", inferred from conformation string," << std::endl
		  << "is not matched by provided sequence:" << std::endl
		  << sequence << ", l = " << sequence.length() << std::endl
		  << "Expected length of " << protein_length << " or "
		  << nuc_length << " (RNA sequence)" << std::endl;
	exit(DATA_ERROR);
    }

    // finally, obtain a codon sequence
    NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());

    // Open translation speed data file and read
    // Make vector of 821, since highest codon is 0x333 = 820
    std::vector<unsigned int> translation_times(821, 0);
    std::vector<int> stop_codon_times;
    ReadTranslationTimes(speedparams_path.c_str(), translation_times.data());

    // Here, we find the three translation times that correspond to
    // the stop codons.
    for (auto it=STOP_CODONS.begin(); it!=STOP_CODONS.end(); ++it)
    {
	stop_codon_times.push_back(translation_times.at(*it));
    }

    // It's time to give the people some information
    if (!g_world_rank)
    {
	*outstream
	    << "# folding_evolution" << std::endl
	    << "# input seq : " << sequence << std::endl
	    << "# n_gens : " << n_gens << std::endl
	    << "# pop size : " << population_size << std::endl
	    << "# temperature : " << temperature << std::endl
	    << "# simulations (no. reevaluations): " << g_world_size <<
	    " (" << reevaluation_size << ")" << std::endl
	    << "# rng seed : " << rng_seed << std::endl
	    << "# translation params : " << speedparams_path << std::endl
	    ;
	cout_header(outstream, aa_sequence, nuc_sequence); 
    }
    
    // Simulation-related variables and parameters
    std::string translation_schedule;
    std::unique_ptr<char []> aa_sequence_str(new char[protein_length+1]);
    double old_fitness = 0.001;
    double fitness;
    double old_folded_fraction;
    double folded_fraction;
    double selection;
    double fixation;
    std::vector<int> temp_sequence;
    std::vector<std::string> latfoldvec_command;
    // For reevaluation of the last accepted generation:
    std::vector<std::string> prev_latfoldvec_command;
    // std::vector<std::string> latmaptraj_command;
    std::string output_dir;
    std::string hdf5_output_file;
    int accept;
    int last_accepted_gen = 0;
    int n_gens_without_mutation = 0;

    // Begin running simulation loop
    for (int gen=0; gen<n_gens; gen++)
    {
	// Compose the simulation parameters
	// For now, the final time is taken to be for stop codon UAA
	translation_schedule = make_translation_schedule(
	    translation_times, codon_sequence, stop_codon_times[0]);
	PrintAASequence(
	    aa_sequence_str.get(), aa_sequence.data(), protein_length);
	latfoldvec_command = compose_latfoldvec_command(
	    latpack_path,
	    aa_sequence_str.get(),
	    folded_conformation,
	    translation_schedule,
	    temperature,
	    latfold_output_frequency);

	if (proc_does_reevaluation && gen)
	{
	    output_dir = h5file_path + "/gen"
		+ std::to_string(last_accepted_gen);
	    hdf5_output_file = output_dir + "/" + "sim"
		+ std::to_string(n_gens_without_mutation + 1) + "_"
		+ std::to_string(g_subcomm_rank) + ".h5";
	    run_latfoldvec(prev_latfoldvec_command, hdf5_output_file);
	    folded_fraction = evaluate_folded_fraction(hdf5_output_file);

	    // Now update old fitness
	    old_folded_fraction = reaverage_folded_fraction(
		old_folded_fraction, folded_fraction, g_subcomm_size,
		g_world_size, n_gens_without_mutation);
	    old_fitness = calculate_fitness(old_folded_fraction);
	}
	else if (!proc_does_reevaluation)
	{
	    output_dir = h5file_path + "/gen"
		+ std::to_string(gen);

	    if (!g_subcomm_rank)
	    {
		if (mkdir(output_dir.c_str(), 0754) == -1)
		{
		    std::cerr << "Error making dir: " << output_dir << std::endl;
		    exit(IO_ERROR);
		}
	    }
	    hdf5_output_file = output_dir + "/" + "sim0_"
		+ std::to_string(g_subcomm_rank) + ".h5";

	    run_latfoldvec(latfoldvec_command, hdf5_output_file);
	    folded_fraction = evaluate_folded_fraction(hdf5_output_file);
	    fitness = calculate_fitness(folded_fraction);
	}
	MPI_Bcast(&fitness, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&folded_fraction, 1, MPI_DOUBLE, g_world_size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&old_fitness, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	selection = (fitness - old_fitness) / old_fitness;
	if (selection == 0.)
	    fixation = 1 / (double)population_size;
	else
	    fixation = (1 - exp(-2 * selection)) /
		(1 - exp(-2 * population_size * selection));

	// Decide whether to accept
	if (!g_world_rank)
	{
	    // Root node decide accept / reject
	    if (fixation > threefryrand())
		accept = true;
	    else
		accept = false;

	    // Output information
	    cout_state(outstream, gen, aa_sequence, nuc_sequence, old_fitness,
		       fitness, accept);
	}

	// Let everyone know accept / reject
	MPI_Bcast(&accept, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	// Update / revert
	if (accept)
	{
	    // Update so that nuc_sequence and fitness are fixed
	    prev_nuc_sequence = nuc_sequence;
	    prev_latfoldvec_command = latfoldvec_command;
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
		temp_sequence = nuc_sequence;
		int m = PointMutateNucSequence(nuc_sequence.data(), nuc_length);
		if (m == mutation_mode)
		    break;
		if (mutation_mode == MutateAll && m >= 0)
		    break;
		nuc_sequence = temp_sequence;
	    } while (true);
	}

	MPI_Bcast(nuc_sequence.data(), nuc_sequence.size(), MPI_INT,
		  0, MPI_COMM_WORLD);
	
	// and update sequences
	NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());
    }
    
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
    const int & final_time)
{
    std::ostringstream os;
    std::vector<std::string> schedule;
    // Note that we skip the first five codons, since they are present
    //   at the start of the MC simulation.
    for (auto codon=codon_sequence.begin()+5;
	 codon!=codon_sequence.end(); ++codon)
    {
	os << times.at(*codon) << " ";
    }
    os << final_time;
    return os.str();
}


// Build a vector of strings that are arguments to invoke
// latFoldVec. Note arguments -seed and -outFile are not added here.
std::vector<std::string> compose_latfoldvec_command(
    const std::string& latpack_path,
    const std::string& aa_sequence,
    const std::string& folded_conformation,
    const std::string& translation_schedule,
    const double& temperature,
    const int& output_frequency)
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
    command.push_back("-ribosomeRelease"); // do we want this?
    // command.push_back("-out=N");
    command.push_back("-out=S");
    command.push_back("-outFreq=" + std::to_string(output_frequency));
    
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
    std::vector<std::string> & latfoldvec_command,
    const std::string & outfile)
{
    // Add additional parameters
    latfoldvec_command.push_back(
	"-seed=" + std::to_string(threefryrand_int() % INT32_MAX));
    latfoldvec_command.push_back("-outFile=" + outfile);
    // latfoldvec_command.push_back("-title=?");

    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
	latfoldvec_command);

    pid_t pid = fork();
    if (pid == -1)
    {
	std::cerr << "Failed to fork" << std::endl;
	exit(IO_ERROR);
    }
    else if (pid == 0)
    {
	// child
	// suppress std out
	int fd = open("/dev/null", O_WRONLY);
	dup2(fd, 1);

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
		exit(LATPACK_ERROR);
	    }
	}
    }
}


// Build a vector of strings that are arguments to invoke
// latMapTraj. Note arguments -traj is not added here.
std::vector<std::string> compose_latmaptraj_command(
    const std::string& latpack_path,
    const std::string& folded_conformation)
{
    std::vector<std::string> command;
    command.push_back(latpack_path + "bin/latMapTraj");
    command.push_back("-ref=" + folded_conformation);
    command.push_back("-lat=CUB");
    command.push_back("-hdf5");
    return command;
}


// FIXME (update to full MPI)
// Use fork and exec to run n_simulations of latMapTraj to analyze the
// latFoldVec simulations that were run.
// this has not been updated yet
void run_latmaptraj(
    std::vector<std::string> & latmaptraj_command,
    const std::string& h5file_base,
    int n_simulations)
{
    for (int i=g_world_rank*n_simulations; i<(g_world_rank+1)*n_simulations; ++i)
    {
	pid_t pid = fork();
	if (pid == -1)
	{
	    std::cerr << "Failed to fork" << std::endl;
	    exit(IO_ERROR);
	}
	else if (pid == 0)
	{
	    // child
	    // we land here for each instance
	    // Add additional parameters
	    latmaptraj_command.push_back("-traj=" + h5file_base +
					 std::to_string(i) + ".h5");
	    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
		latmaptraj_command);

	    // suppress std out
	    int fd = open("/dev/null", O_WRONLY);
	    dup2(fd, 1);

	    execv(cstring_command_vec[0], cstring_command_vec.data());
	}
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
		exit(LATPACK_ERROR);
	    }
	}
    }    
}



// Analyze latFoldVec simulations to determine a fitness value for
// protein under evaluation.
double evaluate_folded_fraction(
    const std::string& filename)
{
    int folded = 0;
    int total_n_folded = 0;

    // Now for each hdf5 file, we need to read in the data
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_id = H5Gopen2(file_id, "traj1", H5P_DEFAULT);
    hid_t attribute_id = H5Aopen(group_id, "Found final struct",  H5P_DEFAULT);

    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 4);

    char found_final[4];

    H5Aread(attribute_id, strtype, found_final);

    if (std::strncmp(found_final, "yes", 3) == 0)
    {
	folded = 1;
    }
	
    H5Tclose(strtype);
    H5Aclose(attribute_id);
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
    
    return folded_fraction;
}


// Calculate fitness from folded fraction.
double calculate_fitness(
    double folded_fraction,
    double f_0)
{
    return 0.001 + folded_fraction / (folded_fraction + f_0);    
}


// Update fitness to incorporate new fitness evaluations.
double reaverage_folded_fraction(
    double old_folded_fraction,
    double new_folded_fraction,
    int n_reevaluators,
    int n_total,
    int n_gens_without_mutation)
{
    double old_factor = (n_total + n_reevaluators * (n_gens_without_mutation - 1))
	/ (double)(n_total + n_gens_without_mutation * n_reevaluators);
    double new_factor = n_reevaluators
	/ (double)(n_total + n_gens_without_mutation * n_reevaluators);
    return old_folded_fraction * old_factor
	+ new_folded_fraction * new_factor;
}


// Print the header for state output
void cout_header(
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
	std::setw(6) << "Accept" << std::endl;

    int line_length = 4;
    line_length += aa_sequence.size() + 1;
    line_length += nuc_sequence.size() + 1;
    line_length += 16 + 2;
    line_length += 6;

    *outstream << "# " << std::string(line_length, '-')<< std::endl;

    *outstream << std::fixed << std::setprecision(2);
    
    return;
}


// Print out what happened each generation.
void cout_state(
    std::ostream * outstream,
    int generation,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double old_fitness,
    double new_fitness,
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
	       << std::setw(6) << ((accept) ? "yes" : "no")
	       << std::endl;

    return;    
}
