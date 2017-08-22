////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_PSSM
#define INCLUDED_PSSM

#include <iostream>
#include <vector>

#include "util.h"

////////////////////////////////////////////////////////////////////////////////
// position-specific-search matrix classes
////////////////////////////////////////////////////////////////////////////////
//// one sequence position of the position-specific search matrix
class PssmPos {

	public:
		PssmPos() : siteindex_(0) {}

		PssmPos( int _index ) : siteindex_( _index ) {}

		void add_weight( float weight );
		std::vector< float > const & weights() const { return weights_; }
//		float weights( unsigned typeindex ) const { return weights_[typeindex]; }
		void scale( float factor );
		float bestweight() const;
		float worstweight() const;
		void print_weights( std::ostream & out ) const;

	private:
		int siteindex_;
		std::vector< float > weights_;

};

std::ostream & operator << ( std::ostream & out, PssmPos const & pssm_pos );

////////////////////////////////////////////////////////////////////////////////
//// the full PS matrix
class PSSM {

	public:
		PSSM() : length_(0), outputlevel_(NORMAL) {}

		void setup( std::string const & filename, bool invert, OutputLevel level = NORMAL );
		void setup( std::string const & target);

		unsigned length() const { return length_; }
		int priority( unsigned siteindex ) const { return priority_[siteindex]; }
		void print( std::ostream & out = std::cout ) const;
		float bestcase( unsigned siteindex ) const { return best_cases_[siteindex]; }
		float score( int siteindex, char letter ) const;
		float bestweight( int siteindex ) const { return positions_[siteindex].bestweight(); }

	private:
		void parse_key( std::string const & line );
		void readfile( std::string const & filename, bool invert );
		void set_priority_and_best_cases();

	private:
		std::vector< int > priority_;
		std::vector< char > key_; // letter to number conversion for weight indices
		std::vector< PssmPos > positions_;
		unsigned length_;
		std::vector< float > best_cases_;
		OutputLevel outputlevel_;

};

std::ostream & operator << ( std::ostream & out, PSSM const & pssm );

#endif
