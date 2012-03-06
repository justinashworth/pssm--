////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_Sequence
#define INCLUDED_Sequence

#include <iostream>
#include <vector>

#include "util.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
class Gene {

	public:
		//// constructors
		Gene()
			: finalized_(false)
		{}

		Gene( char hdr[] ) { header_.assign( hdr ); finalized_ = false; }

		std::vector< char > const & sequence() const { return sequence_; }
		std::string const & header() const { return header_; }
		void readline( char linebuff[], int const chars );

		// iterators to provide read access to the gene sequence
		std::vector< char >::const_iterator begin() const { return sequence_.begin(); }
		std::vector< char >::const_iterator end() const { return sequence_.end(); }

		void finalize();

		//// abbreviated summary of the gene sequence
		void print( std::ostream & out = std::cout ) const;

	private:
		// try to make this data private...
		std::vector< char > sequence_; // the meat
		std::string header_;
		bool finalized_;
};

// output stream operator for Gene
std::ostream & operator<< ( std::ostream & out, Gene const & gene );

////////////////////////////////////////////////////////////////////////////////////////////////////
/// container and manager for a list of genes
class GeneList {

	public:
		GeneList()
			: numseqs_(0),
				numbps_(0),
				outputlevel_(NORMAL)
		{}

		GeneList( std::string const filename, OutputLevel level = NORMAL )
			: numseqs_(0),
				numbps_(0),
				outputlevel_(level)
		{
			readfile( filename );
		}

		unsigned numseqs() const { return numseqs_; }
		unsigned numbps() const { return numbps_; }
		// iterators to provide read access to the gene sequence list
		std::vector< Gene >::const_iterator begin() const { return genes_.begin(); }
		std::vector< Gene >::const_iterator end() const { return genes_.end(); }
		void print( std::ostream & out = std::cout ) const;

	private:
		void finalize();

		//// opens the input file and feeds lines to the Gene class
		void readfile( std::string const filename );

	private:
		std::vector< Gene > genes_;
		unsigned numseqs_, numbps_;
		OutputLevel outputlevel_;
};

std::ostream & operator << ( std::ostream & out, GeneList const & genelist );

#endif
