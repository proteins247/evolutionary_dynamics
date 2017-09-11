/* 
 * folding_evolution
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <csignal> 		// eventually need to think about handling interrupts
#include <cmath>
#include <map>

#include <getopt.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#define CPLUSPLUS __cplusplus
#undef __cplusplus
// When compiling Random123 header files for rng.h, we don't want c++
// features (they're not compatible with extern "C"), so we
// temporarily undefine __cplusplus, saving its value in CPLUSPLUS.
extern "C" {
#include "../gencode.h"
#include "../rng.h"
}
#define __cplusplus CPLUSPLUS
#undef CPLUSPLUS


static const std::string helptext =
    "folding_evolution\n"
    "Usage: folding_evolution [OPTIONS] SEQUENCE N_GENS \\\n"
    "         --native-fold NFOLD --speed-params FILE\n"
    "Carry out folding evolutionary simulations\n"
    "\n"
    "  -n, --native-fold=CONF    Required. The native conformation of the\n"
    "                            as a latPack-style move sequence.\n"
    "  -s, --speed-params=FILE   Required. Path to file with translation\n"
    "                            times of the codons.\n"
    "  -d, --tempfile-path=PATH  LatFoldVec trajectory files written to this\n"
    "                            directory. Default=./\n"
    "  -l, --latpack-path=PATH   Latpack binaries here. Default=\n"
    "                            /n/home00/vzhao/pkg/latPack/1.9.1-6/\n"
    "  -N, --sims-per-gen=N      Run N instances of latFoldVec per\n"
    "                            generation. Need N processors available.\n"
    "  -o, --outfile=FILE        Output will be sent to this file instead of\n"
    "                            STDOUT.\n"
    "  -p, --population-size=N   Size of population for evolutionary\n"
    "                            dynamics. Default=100\n"
    "      --random-codons       If SEQUENCE is of amino acids, corresponding\n"
    "                            RNA sequence codons will be randomly picked.\n"
    "  -r, --seed=N              RNG seed. Default=1\n"
    "  -t, --temperature=T       Temperature of latFoldVec simulations.\n"
    "                            Default=0.3\n"
    "      --help   Display this help and exit.\n"
    ;

// error values
static const int PARSE_ERROR = 1;
static const int DATA_ERROR = 2;
static const int IO_ERROR = 3;
static const int LATPACK_ERROR = 4;

// default values
static const std::string DEFAULT_LATPACK_PATH = "/n/home00/vzhao/pkg/latPack/1.9.1-6/";
static const std::string DEFAULT_TEMPFILE_PATH = "./";
static const uint64_t DEFAULT_SEED = 1;
static const int DEFAULT_N_SIMS = 1;
static const int DEFAULT_LATFOLD_OUTFREQ = 100;
static const int DEFAULT_POPULATION_SIZE = 100;
static const double DEFAULT_TEMPERATURE = 0.3;

typedef std::map<Codon, int> CodonIntMap;


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
    const CodonIntMap & times,
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


// Use fork and exec to run n_simulations of latFoldVec co-translational
// folding simulations.
//
// @param latfoldvec_command The vector of strings built
//        by compose_latfoldvec_command.
// @param h5file_base latFoldVec program output will go to
//        a file whose pathname begins with this.
// @param n_simulations Run this many folding simulations.
void run_latfoldvec(
    std::vector<std::string> & latfoldvec_command,
    const std::string & h5file_base,
    int n_simulations=1);


// Build a vector of strings that are arguments to invoke
// latMapTraj. Note arguments -traj is not added here.
std::vector<std::string> compose_latmaptraj_command(
    const std::string& latpack_path,
    const std::string& folded_conformation);


// Use fork and exec to run n_simulations of latMapTraj to analyze the
// latFoldVec simulations that were run.
void run_latmaptraj(
    std::vector<std::string> & latmaptraj_command,
    const std::string& h5file_base,
    int n_simulations=1);


// Analyze latFoldVec simulations to determine a fitness value for
// protein under evaluation.
double evaluate_fitness(
    const std::string& folded_conformation,
    const std::string& h5file_base,
    int n_simulations=1);


// Print the header for state output
// we would want
void cout_header(
    std::ostream * outstream,
    const std::vector<AminoAcid> & aa_sequence,
    const std::vector<int> & nuc_sequence);


// Print out what happened each generation.
void cout_state(
    std::ostream * outstream,
    int generation,
    std::vector<AminoAcid> & aa_sequence,
    std::vector<int> & nuc_sequence,
    double fitness);


int main(int argc, char** argv)
{
    // Variables to be determined by commandline options: 
    std::string sequence;
    std::string folded_conformation;
    std::string speedparams_path;
    std::string outfile_path;
    std::string latpack_path = DEFAULT_LATPACK_PATH;
    std::string tempfile_path = DEFAULT_TEMPFILE_PATH;
    std::ostream * outstream = &std::cout;
    std::ofstream outfile;
    uint64_t rng_seed = DEFAULT_SEED;

    // how many folding simulations run per generation
    // is controlled by this parameter
    int n_simulations = DEFAULT_N_SIMS;

    // This particular param currently cannot be changed by
    // commandline arguments
    int latfold_output_frequency = DEFAULT_LATFOLD_OUTFREQ;

    int n_gens;			// How many generations to run
    int random_codons = 0;
    int population_size = DEFAULT_POPULATION_SIZE;
    double temperature = DEFAULT_TEMPERATURE;
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
	{"tempfile-path", required_argument, NULL, 'd'},
	{"latpack-path", required_argument, NULL, 'l'},
	{"native-fold", required_argument, NULL, 'n'},
	{"sims-per-gen", required_argument, NULL, 'N'},
	{"random-codons", no_argument, &random_codons, 1},
	{"outfile", required_argument, NULL, 'o'},
	{"population-size", required_argument, NULL, 'p'},
	{"seed", required_argument, NULL, 'r'},
	{"speed-params", required_argument, NULL, 's'},
	{"temperature", required_argument, NULL, 't'},
	{"help", no_argument, NULL, 'h'},
	{NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "d:l:n:N:o:p:r:s:t:h",
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
	    tempfile_path = optarg;
	    break;
	case 'l':
	    latpack_path = optarg;
	    break;
	case 'n':
	    folded_conformation = optarg;
	    break;
	case 'N':
	    try
	    {
		n_simulations = std::stoi(optarg);
	    }
	    catch (...)
	    {
		std::cerr << "Failed to convert n_simulations: " << optarg
			  << std::endl;
		exit(PARSE_ERROR);
	    }
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
		std::cerr << "Failed to convert ppoulation size: " << optarg
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

    // initialize RNG
    set_threefry_array(rng_seed);

    // Setup output
    if (!outfile_path.empty())
    {
	// You can think of the contents in this block as placeholder
	// in case one day there is a switch to non-text based output
	// because the user can just as easily redirect stdout

	outfile.open(outfile_path);
	if (!outfile)
	{
	    std::cerr << "Problem with output file: "
		      << outfile_path << std::endl;
	    exit(IO_ERROR);
	}
	outstream = &outfile;
    }
    
    // Begin simulation setup
    // Process the user-provided sequence
    std::vector<int> nuc_sequence;
    std::vector<AminoAcid> aa_sequence;
    std::vector<Codon> codon_sequence;

    // Infer lengths from folded conformation length
    unsigned int protein_length = folded_conformation.length() + 1;
    unsigned int nuc_length = protein_length * 3;

    // After inferring length, reserve space
    nuc_sequence.resize(nuc_length);
    aa_sequence.resize(protein_length);
    codon_sequence.resize(protein_length);

    // We use reserve because we directly access the underlying arrays of
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
	AASeqToNucSeq(aa_sequence.data(), nuc_sequence.data(), protein_length,
	    random_codons);
	
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
    std::vector<unsigned int> translation_times_vec(821, 0);
    ReadTranslationTimes(speedparams_path.c_str(), translation_times_vec.data());

    // This is hackish
    CodonIntMap translation_times;
    int idx = 0;
    for (auto it=translation_times_vec.begin(); it!=translation_times_vec.end();
	 ++it, ++idx)
    {
	if (*it != 0)
	    translation_times.insert(std::pair<Codon, int>((Codon)idx, *it));
    }

    // It's time to give the people some information
    *outstream
	<< "# folding_evolution" << std::endl
	<< "# input seq : " << sequence << std::endl
	<< "# n_gens : " << n_gens << std::endl
	<< "# pop size : " << population_size << std::endl
	<< "# temperature : " << temperature << std::endl
	<< "# sims per gen : " << n_simulations << std::endl
	<< "# rng seed : " << rng_seed << std::endl
	;

    cout_header(outstream, aa_sequence, nuc_sequence); 
    
    // Simulation-related variables and parameters
    std::string translation_schedule;
    std::unique_ptr<char []> aa_sequence_str(new char[protein_length+1]);
    double old_fitness = 0.001;
    double fitness;
    double selection;
    double fixation;
    std::vector<int> old_nuc_sequence = nuc_sequence;
    std::vector<std::string> latfoldvec_command;
    std::vector<std::string> latmaptraj_command;
    std::string outfile_base;
    std::string latpack_outfile_signature = std::to_string(
	threefryrand_int() % 1000);

    // Begin running simulation
    for (int gen=0; gen<n_gens; gen++)
    {
	// Compose the simulation parameters
	translation_schedule = make_translation_schedule(
	    translation_times, codon_sequence);

	PrintAASequence(aa_sequence_str.get(), aa_sequence.data(), protein_length);

	latfoldvec_command = compose_latfoldvec_command(
	    latpack_path,
	    aa_sequence_str.get(),
	    folded_conformation,
	    translation_schedule,
	    temperature,
	    latfold_output_frequency);

	outfile_base = tempfile_path + "/fold_evo_sim__gen"
	    + std::to_string(gen) + "_"
	    + latpack_outfile_signature + "_sim";

	// run simulation(s)
	run_latfoldvec(latfoldvec_command, outfile_base, n_simulations);

	// run analysis
	latmaptraj_command = compose_latmaptraj_command(
	    latpack_path, folded_conformation);

	run_latmaptraj(latmaptraj_command, outfile_base, n_simulations);

	// evaluate fitness
	fitness = evaluate_fitness(
	    folded_conformation,
	    outfile_base,
	    n_simulations);
	selection = (fitness - old_fitness) / old_fitness;
	fixation = (1 - exp(-2 * selection)) /
	    (1 - exp(-2 * population_size * selection));

	// decide whether to accept
	if (fixation > threefryrand())
	{
	    old_nuc_sequence = nuc_sequence;
	}
	else
	{
	    nuc_sequence = old_nuc_sequence;
	}

	// Output information
	cout_state(outstream, gen, aa_sequence, nuc_sequence, fitness);

	// Make new mutation
	while (PointMutateNucSequence(nuc_sequence.data(), nuc_length) < 0);

	// and update sequences
	NucSeqToCodonSeq(nuc_sequence.data(), nuc_length, codon_sequence.data());
	NucSeqToAASeq(nuc_sequence.data(), nuc_length, aa_sequence.data());

    }

    return 0;
} // end main


// Print this program's help message
void print_help()
{
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
    const CodonIntMap & times,
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
    command.push_back("-out=S");
    command.push_back("-outFreq=" + std::to_string(output_frequency));
    
    // Previously we did the what is commented out, but we want to
    // keep the command as a vector.
    // std::ostringstream os;
    // std::copy(command.begin(), command.end(),
    // std::ostream_iterator<std::string>(os, " ")); return os.str();
    return command;
}


// Use fork and exec to run n_simulations of latFoldVec co-translational
// folding simulations.
void run_latfoldvec(
    std::vector<std::string> & latfoldvec_command,
    const std::string & h5file_base,
    int n_simulations)
{
    for (int i=0; i<n_simulations; ++i)
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
	    latfoldvec_command.push_back(
		"-seed=" + std::to_string(threefryrand_int() % INT32_MAX));
	    latfoldvec_command.push_back(
		"-outFile=" + h5file_base + std::to_string(i) + ".h5");
	    // latfoldvec_command.push_back("-title=?");

	    std::vector<char *> cstring_command_vec = string_vec_to_cstring_vec(
		latfoldvec_command);

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


// Use fork and exec to run n_simulations of latMapTraj to analyze the
// latFoldVec simulations that were run.
void run_latmaptraj(
    std::vector<std::string> & latmaptraj_command,
    const std::string& h5file_base,
    int n_simulations)
{
    for (int i=0; i<n_simulations; ++i)
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
double evaluate_fitness(
    const std::string& folded_conformation,
    const std::string& h5file_base,
    int n_simulations)
{
    // Now for each hdf5 file, we need to read in the data

    // temporary placeholder, just return a random number lol
    return threefryrand();
}


// Print the header for state output
void cout_header(
    std::ostream * outstream,
    const std::vector<AminoAcid> & aa_sequence,
    const std::vector<int> & nuc_sequence)
{
    auto shorten = [](std::string str, int len)
    {
	return str.substr(0, len);
    };

    *outstream << std::endl;
    *outstream << std::left <<
	"# Gen|" <<
	std::setw(aa_sequence.size()) <<
	shorten("AA sequence", aa_sequence.size()) << "|" <<
	std::setw(nuc_sequence.size()) <<
	shorten("Nuc sequence", nuc_sequence.size()) << "|" <<
	std::setw(5) << "Fit." <<
	std::endl;

    int line_length = 6;
    line_length += aa_sequence.size() + 1;
    line_length += nuc_sequence.size() + 1;
    line_length += 5;

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
    double fitness)
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
	       << std::setw(5) << fitness
	       << std::endl;

    return;    
}


