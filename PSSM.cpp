////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <math.h> // fabs
#include <sstream>
#include <vector>
#include <cstdlib> // exit, EXIT_FAILURE
#include <algorithm> // std::sort

#include "PSSM.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// position-specific-search matrix classes
////////////////////////////////////////////////////////////////////////////////
//// one sequence position of the position-specific search matrix
void
PssmPos::add_weight( float weight )
{
	weights_.push_back( weight );
}

void
PssmPos::scale( float factor )
{
	for ( std::vector< float >::iterator weight( weights_.begin() );
	      weight != weights_.end(); ++weight ) { *weight *= factor; }
}

float
PssmPos::bestweight() const
{
	float best( weights_[0] );
	for ( unsigned i(1); i < weights_.size(); ++i ) {
		if ( weights_[i] < best ) best = weights_[i];
	}
	return best;
}

float
PssmPos::worstweight() const
{
	float worst( weights_[0] );
	for ( unsigned i(1); i < weights_.size(); ++i ) {
		if ( weights_[i] > worst ) worst = weights_[i];
	}
	return worst;
}

void
PssmPos::print_weights( std::ostream & out ) const
{
	out << siteindex_;
	for ( std::vector< float >::const_iterator weight( weights_.begin() );
	      weight != weights_.end(); ++weight ) { out << " " << *weight; }
}

std::ostream & operator << ( std::ostream & out, PssmPos const & pssm_pos )
{
	pssm_pos.print_weights( out );
	return out;
}

////////////////////////////////////////////////////////////////////////////////
//// the full PS matrix

void
PSSM::setup(
	std::string const & filename,
	bool invert,
	OutputLevel level // = NORMAL
)
{
	outputlevel_ = level;
	readfile( filename, invert );
}

void
PSSM::setup(std::string const & target)
{
	int siteindex(0);
	std::list<char> const nucs( nucleotides() );

	for(std::list<char>::const_iterator nt(nucs.begin()); nt!=nucs.end(); ++nt){
		key_.push_back(*nt);
	}

	for(std::string::const_iterator it(target.begin()); it!=target.end(); ++it){
		PssmPos new_pssm_pos(siteindex);
		// all of the nucs specified by code letter
		std::list<char> const deg( degen(*it) );
		// check each of the key_ nucs for representation in deg
		for(std::vector<char>::const_iterator k(key_.begin()); k!=key_.end(); ++k){
			float weight(0);
			for(std::list<char>::const_iterator b(deg.begin()); b!=deg.end(); ++b){
				if (*k == *b) weight = -1;
			}
			new_pssm_pos.add_weight(weight);
		}
		positions_.push_back( new_pssm_pos );
	}
	length_ = positions_.size();

	set_priority_and_best_cases();
}

void
PSSM::print(
	std::ostream & out
) const
{
	out << "key";
	for ( std::vector< char >::const_iterator letter( key_.begin() );
	      letter != key_.end(); ++letter ) { out << " " << *letter; }
	out << std::endl;

	for ( std::vector< PssmPos >::const_iterator pssm_pos( positions_.begin() );
  pssm_pos != positions_.end(); ++pssm_pos ) { out << *pssm_pos << std::endl; }
}

float
PSSM::score(
	int siteindex,
	char letter
) const
{
//	std::cout << "DEBUG: letter=" << letter << std::endl;
	for ( unsigned i(0); i < key_.size(); ++i ) {
		if ( letter == key_[i] ) {
			return positions_[ siteindex ].weights()[i];
		}
	}
	std::cerr << "ERROR: score failed for unknown character (" << letter << ") at " << siteindex << std::endl;
//	exit(EXIT_FAILURE);
	return 0;
}

void
PSSM::parse_key( std::string const & line )
{
	std::istringstream linestream( line );
	std::string dummy;
	linestream >> dummy; // "key"
	char letter;
	while ( linestream >> letter ) key_.push_back( upper( letter ) );
	if(outputlevel_ > NORMAL){
		std::cout << "key:" << std::endl;
		for(std::vector<char>::const_iterator k(key_.begin()); k!=key_.end(); ++k){
			std::cout << *k << std::endl;
		}
	}
}

void
PSSM::readfile(
	std::string const & filename,
	bool invert
)
{
	// open input file
	std::ifstream file;
	file.open( filename.c_str() );
	if ( !file ) {
		std::cerr << "ERROR: unable to open PSSM file " << filename.c_str() << std::endl;
		exit(EXIT_FAILURE);
	}
	if ( outputlevel_ >= NORMAL ) std::cout << "Reading PSSM file " << filename << std::endl;

	// read input file
	std::string line;
	while ( getline( file, line ) ) {
		if ( line[0] == '#' ) continue;
		if ( line.substr(0,3) == "key" || line.substr(0,3) == "KEY" ) parse_key( line );
		else {
			std::istringstream linestream( line );
			int siteindex;
			linestream >> siteindex;
			PssmPos new_pssm_pos(siteindex);
			float weight;
			while ( linestream >> weight ) {
				if ( new_pssm_pos.weights().size() == key_.size() ) {
					std::cerr << "Error: more weights given for " << siteindex
					          << " than denoted in key!" << std::endl;
					exit(EXIT_FAILURE);
				}
				new_pssm_pos.add_weight( weight );
			}
			if ( new_pssm_pos.weights().size() < key_.size() ) {
				std::cerr << "Error: less weights given for " << siteindex
				          << " than denoted in key!" << std::endl;
				exit(EXIT_FAILURE);
			}
			positions_.push_back( new_pssm_pos );
		}
	}
	length_ = positions_.size();

	if ( invert ) {
		for ( std::vector< PssmPos >::iterator pssm_pos( positions_.begin() );
		      pssm_pos != positions_.end(); ++pssm_pos ) { pssm_pos->scale( -1 ); }
	}
	if ( outputlevel_ >= MINIMAL ) print();
	set_priority_and_best_cases();
}

//// this sets up fairly important optimizations of the naive approach
void
PSSM::set_priority_and_best_cases()
{
	unsigned numpos( positions_.size() );

	// prioritizing the scoring of the positions with the highest variability, because these are the positions with the highest information content
	// compile list of ranges
	std::vector< std::pair< unsigned, float > > ranges;
	for ( unsigned i(0); i < numpos; ++i ) {

		float range( fabs( positions_[i].bestweight() -
		                   positions_[i].worstweight() ) );

		ranges.push_back( std::pair< unsigned, float >( i, range ) );
	}
	// sorted list of vector indices by range (descending)
	std::sort( ranges.begin(), ranges.end(), secondfloatdesc );

	// priority to be a vector of weight vector indices ordered by importance
	priority_.assign( numpos, 0 ); // predimension
	for ( unsigned i(0); i < numpos; ++i ) priority_[i] = ranges[i].first;

	// setup for best-case-scenario early rejection
	// for a given position, this gives the best possible score contributed by the remaining positions
	// indices must correspond to priority vector
	best_cases_.assign( numpos, 0. ); // predimension
	for ( unsigned i(0); i < numpos; ++i ) {
		for ( unsigned j(i+1); j < numpos; ++j ) {
			best_cases_[i] += positions_[ priority_[j] ].bestweight();
		}
	}

	// debug output for priority and best cases
	// probably want to implement a more general debug output scheme for PSSM class with own method
	if ( outputlevel_ >= VERBOSE ) {
		std::cout << "\nPSSM debug output:\nOptimal priorities for the order of scoring positions:\n";
		for ( unsigned i(0), size( priority_.size() ); i < size; ++i ) {
			std::cout << i << " " << priority_[i] << " " << ranges[i].second << '\n';
		}
		std::cout << "\nBest additional score possible given remaining positions:\n";
		for ( unsigned i(0), size( best_cases_.size() ); i < size; ++i ) {
			std::cout << i << " " << best_cases_[i] << '\n';
		}
		std::cout << std::endl;
	}

}

std::ostream & operator << ( std::ostream & out, PSSM const & pssm )
{
	pssm.print( out );
	return out;
}

