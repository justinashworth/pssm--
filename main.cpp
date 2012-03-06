////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////
//
// c++ version of the position-specific-matrix sequence search.
// MUCH faster, albeit slower to develop then the Perl or Python versions.
//

/* Status:

2008-12
	// cleaned up classes, split them into separate files w/ headers, added options
	// searching for top fifty 24 bp sites in typical whole genome sequence takes a few minutes

2007-07-24
	// compiled with -O3 and without -ggdb symbols
	0.116u 0.001s 0:00.12 91.6%     0+0k 0+0io 0pf+0w

2007-07-05
	// with fwd and rvs scoring, early rejection, scoring priority.
	0.931u 0.003s 0:00.93 100.0%    0+0k 0+0io 0pf+0w

	// with fwd and rvs scoring, and early rejection
	4.467u 0.005s 0:04.47 99.7%     0+0k 0+0io 0pf+0w

	// with fwd and rvs scoring
	8.902u 0.006s 0:08.91 99.8%     0+0k 0+0io 0pf+0w

2007-07-04
	Basic functionality implemented.

*/

////////////////////////////////////////////////////////////////////////////////

#include <dirent.h>
#include <fstream>
#include <iostream>
#include <list>
#include <cstdlib> // exit, EXIT_FAILURE

#include "util.h"
#include "Sequence.h"
#include "Hits.h"
#include "PSSM.h"
#include "TargetSearch.h"

////////////////////////////////////////////////////////////////////////////////
void usage_error()
{
	std::cerr << "\n"
	 << " -s|--seq|--sequence     sequencefile   : FASTA format\n"
	 << " -l|--list               seqlistfile    : file with list of FASTA files\n"
	 << " -p|--pssm               pssmfilename   : weight matrix\n"
	 << " -inv                                   : invert weights (for positive weights)\n"
	 << " -n|--numhits|--hits     #              : number of hits (20)\n"
	 << " -v|--verbose                           : more output\n"
	 << " -m|--minimal|--mute                    : less output\n"
	 << "example: [executable] -s genes.dna -p mso-xray.pssm\n"
	 << "\n";
	exit(EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
int main( int argc, char *argv[] ) {

	std::cout << std::endl;

	std::string seqfilename, seqlistname, pssmfilename;
	unsigned numhits(20);
	bool invert_pssm(false);
	OutputLevel outputlevel(NORMAL);

	// parse command line arguments
	for ( int i(1); i < argc; ++i ) {

		std::string arg( argv[i] );

		if ( arg == "-s" || arg == "--seq" || arg == "--sequence" ) {
			if ( ++i >= argc ) usage_error();
			seqfilename = argv[i];

		} else if ( arg == "-l" || arg == "--list" ) {
			if ( ++i >= argc ) usage_error();
			seqlistname = argv[i];

		} else if ( arg == "-p" || arg == "--pssm" ) {
			if ( ++i >= argc ) usage_error();
			pssmfilename = argv[i];

		} else if ( arg == "-n" || arg == "--numhits" || arg == "--hits" ) {
			if ( ++i >= argc ) usage_error();
			numhits = atoi( argv[i] );

		} else if ( arg == "-inv" ) {
			invert_pssm = true;

		} else if ( arg == "-v" || arg == "--verbose" ) {
			outputlevel = VERBOSE;

		} else if ( arg == "-m" || arg == "--minimal" || arg == "--mute" ) {
			outputlevel = MINIMAL;

		} else if ( arg == "-h" || arg == "--help" ) {
			usage_error();

		}
	}

	// get sequence filenames
	std::list< std::string > filenames;

	// single filename given
	if ( !seqfilename.empty() ) filenames.push_back( seqfilename );

	// file with list of filenames given
	else if ( !seqlistname.empty() ) {
		std::ifstream seqlistfile;
		seqlistfile.open( seqlistname.c_str() );
		if ( !seqlistfile ) {
			std::cerr << "ERROR: couldn't open " << seqlistname << std::endl;
			exit(EXIT_FAILURE);
		}
		std::string line;
		while ( getline( seqlistfile, line ) ) {
			filenames.push_back( line );
		}

	// if all else fails, try to find and use sequence files in current directory
	} else {
		DIR *dp( opendir(".") );
		struct dirent *dirp;
		while ( ( dirp = readdir(dp) ) ) {
			std::string name( dirp->d_name );
			unsigned const l( name.size() );
			// check file extension
			unsigned i(l-1);
			while ( i > 0 ) {
				char ichar( name.at(--i) );
				if ( ichar == '.' ) break;
			}
			if ( i == 0 ) continue;
			std::string ext( name.substr(i) );
			if ( ext != ".fa" && ext != ".dna" && ext != ".seq" ) continue;
			filenames.push_back( std::string( dirp->d_name ) );
		}
		closedir(dp);
	}

	if ( pssmfilename.empty() ) {
		std::cout << "ERROR: No pssm file given." << std::endl;
		usage_error();
	}
	TargetSearch search( pssmfilename, numhits, invert_pssm, outputlevel );
	// perform the search, operates as a functor over gene files
	for ( std::list< std::string >::const_iterator name( filenames.begin() );
	      name != filenames.end(); ++name ) {
		search.scan_seq( *name );
	}
	search.print_results();

	return 0;
}

