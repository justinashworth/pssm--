////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <algorithm> // std::transform

#include "Sequence.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

//// read and filter a line of input sequence from the file stream
void Gene::readline( char linebuff[], int const chars )
{
	// use temporary filter container to avoid constly multiple redimensioning
	std::vector< char > clean( chars );
	int nclean(0);
	// filter as rapidly as possible
	for ( int i(0); i < chars; ++i ) {
		char const thischar = linebuff[i];
		if ( isnuc(thischar) ) clean[nclean++] = thischar;
		else std::cerr << header_ << ": unrecognized letter (" << thischar << ") at position " << i << std::endl;
	}
	sequence_.insert( sequence_.end(), clean.begin(), clean.begin() + nclean );
}

void
Gene::finalize()
{
	if ( finalized_ ) return;
	// ensure all upper-case
	std::transform( sequence_.begin(), sequence_.end(), sequence_.begin(), upper );
	finalized_ = true;
}

//// abbreviated summary of the gene sequence
void
Gene::print( std::ostream & out ) const
{
	unsigned maxhdr(15), maxseq(40);
	if ( header_.size() < maxhdr ) out << header_;
	else out << header_.substr(0,15);
	out << ": ";
	unsigned seqsize( sequence_.size() );
	if ( seqsize <= maxseq ) out << sequence_;
	else {
		for ( unsigned i(0); i < maxseq/2; ++i ) out << sequence_[i];
		out << " ... ";
		for ( unsigned i(maxseq/2); i > 0; --i ) out << sequence_[seqsize-i];
	}
	out << " (" << seqsize << " bp)" << std::endl;
}

// output stream operator for Gene
std::ostream & operator<< ( std::ostream & out, Gene const & gene )
{
	gene.print( out );
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// container and manager for a list of genes
void GeneList::finalize()
{
	if ( genes_.size() != 0 ) {
		for ( std::vector< Gene >::iterator gene( genes_.begin() );
		      gene != genes_.end(); ++gene ) {
			gene->finalize();
		}
	}
	else std::cerr << "ERROR: no genes in the list!" << std::endl;

	numseqs_ = genes_.size();
	numbps_ = 0;
	for ( std::vector< Gene >::const_iterator gene( genes_.begin() );
  gene != genes_.end(); ++gene ) {
		numbps_ += gene->sequence().size();
	}
}

//// opens the input file and feeds lines to the Gene class
void GeneList::readfile( std::string const filename )
{
	std::ifstream file;
	file.open( filename.c_str() );
	if ( outputlevel_ >= NORMAL ) std::cout << "\nReading sequence file " << filename << std::endl;

	// read input file
	// using a char array is faster than reading in/parsing strings
	char linebuff[4096];
	bool header_read(false);
	while ( file.getline( linebuff, 4096 ) ) {
//		std::cout << "DEBUG: " << linebuff << std::endl;

		if ( linebuff[0] == '>' ) { // FASTA header
//			std::cout << "DEBUG: fasta label: " << linebuff << std::endl;
			header_read = true;
			genes_.push_back( Gene( linebuff ) );
			continue;
		}
		if ( header_read ) {
			genes_.back().readline( linebuff, file.gcount() );
//			std::cout << "DEBUG: fasta sequence linebuff " << linebuff << std::endl;
		}
	}
	// get the final remaining contents of linebuff
	if ( header_read ) {
		genes_.back().readline( linebuff, file.gcount() );
//		std::cout << "DEBUG: fasta sequence final linebuff " << linebuff << std::endl;
	}
	finalize();
}

void
GeneList::print( std::ostream & out ) const
{
	for ( std::vector< Gene >::const_iterator gene( genes_.begin() );
  gene != genes_.end(); ++gene ) {
		out << *gene;
	}
}

std::ostream & operator << ( std::ostream & out, GeneList const & genelist )
{
	genelist.print( out );
	return out;
}

