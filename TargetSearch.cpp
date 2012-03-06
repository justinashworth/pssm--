////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <vector>

#include "util.h"
#include "TargetSearch.h"

TargetSearch::TargetSearch(
	std::string const & pssmfilename,
	unsigned maxhits,
	bool invert_pssm,
	OutputLevel outputlevel // = NORMAL
)
	: numseqs_(0),
		numbps_(0),
		outputlevel_(outputlevel)
{
	pssm_.setup( pssmfilename.c_str(), invert_pssm, outputlevel );
	hits_.maxhits( maxhits );
	hits_.outputlevel( outputlevel );
}

void
TargetSearch::scan_seq( std::string const & filename )
{
	GeneList genelist( filename, outputlevel_ );
	numseqs_ += genelist.numseqs();
	numbps_ += genelist.numbps();
	for ( std::vector< Gene >::const_iterator gene( genelist.begin() );
	      gene != genelist.end(); ++gene ) {
		scan_seq( *gene );
	}
}

void
TargetSearch::print_results( std::ostream & out ) const
{
	std::cout << std::endl;
	std::cout << numseqs_ << " sequences with a total of "
	          << numbps_ << " basepairs searched." << std::endl;
//	hits_.print( out ); // basic output of hits with no markup

	// more informative output of hits by postponed (re)evaluation
	std::cout << std::showpoint << std::fixed << std::setprecision(2);
	for ( std::list< Hit >::const_iterator h( hits_.hits().begin() ), end( hits_.hits().end() );
	      h != end; ++h ) {
		std::cout << h->score() << " ";
		for ( unsigned i(0), size( h->sequence().size() ); i < size; ++i ) {
			char bp( h->sequence()[i] );
			// the basepair letter is made lowercase if it does not represent the best case
			if ( pssm_.score(i,bp) > pssm_.bestweight(i) ) bp = lower( bp );
			std::cout << bp;
		}
		std::cout << " " << h->source() << " " << h->seqindex();
		if ( h->rvs() ) out << " (rvs)";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void TargetSearch::scan_seq( Gene const & gene )
{
	if ( outputlevel_ >= NORMAL ) {
		std::cout << "Searching gene ";
		gene.print();
	}

	unsigned const dotfreq( 100000 );
	if ( outputlevel_ >= VERBOSE ) {
		std::cerr << "(Each dot represents " << dotfreq << " basepairs searched.)" << std::endl;
	}

	std::vector< char > const & sequence( gene.sequence() );
	unsigned start(0), length( pssm_.length() );
	unsigned end( sequence.size() - length );

	while ( start <= end ) {
		float score(0.), worst( hits_.worst() );

//		std::vector<char> dbg( gene.begin()+start, gene.begin()+start+length );
//		std::cout << dbg << std::endl;

		// score forward site
		for ( unsigned p(0); p < length; ++p ) {
			unsigned i( pssm_.priority(p) ); // pssm position index
			score += pssm_.score( i, sequence[start+i] );
			// early rejection
			if ( hits_.full() && ( score + pssm_.bestcase(p) > worst ) ) break;
			if ( p == length-1 ) { // last position
				if ( hits_.full() && score >= worst ) break;
				std::vector<char> hitseq( gene.begin()+start, gene.begin()+start+length );
				hits_.add_hit( score, hitseq, gene.header(), start );
			}
		}

		// score reverse complement
		// separate from above to take advantage of independent early rejection
		score = 0.; worst = hits_.worst();
		for ( unsigned p(0); p < length; ++p ) { // priority index
			unsigned i( pssm_.priority(p) ); // pssm position index
			// careful here: give FWD index, but reverse complement base
			score += pssm_.score( i, comp( sequence[ start+length-i-1 ] ) );
			if ( hits_.full() && ( score + pssm_.bestcase(p) > worst ) ) break;
			if ( p == length-1 ) {
				if ( hits_.full() && score >= worst ) break;
				std::vector<char> hitseq( gene.begin()+start, gene.begin()+start+length );
				hits_.add_hit( score, hitseq, gene.header(), start, true );
			}
		}

		if ( outputlevel_ >= VERBOSE ) {
			if ( start % dotfreq == 0 ) std::cerr << ".";
		}

		++start;
	}
}

