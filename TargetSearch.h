////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_TargetSearch
#define INCLUDED_TargetSearch

#include <iostream>

#include "Sequence.h"
#include "Hits.h"
#include "PSSM.h"

// the highest-level (application) class
class TargetSearch {

	public:
		TargetSearch(
			std::string const & pssmfilename,
			unsigned maxhits,
			bool invert_pssm,
			OutputLevel outputlevel = NORMAL
		);

		void scan_seq( std::string const & filename );
		void print_results( std::ostream & out = std::cout ) const;

	private: // methods
		void scan_seq( Gene const & gene );

	private: // data
		HitManager hits_;
		PSSM pssm_;
		unsigned numseqs_, numbps_;
		OutputLevel outputlevel_;
};

#endif
