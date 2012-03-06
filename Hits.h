////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_Hits
#define INCLUDED_Hits

#include <iostream>
#include <list>
#include <vector>

#include "util.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
class Hit {

	public:
		Hit()
			: score_( 0.0 ),
				source_( "" ),
				seqindex_( 0 ),
				rvs_( false )
		{}

		Hit(
			std::vector< char > const & _sequence,
			float _score,
			std::string const & _source,
			unsigned _seqindex,
			bool _rvs = false
		)
			: sequence_( _sequence ),
				score_( _score ),
				source_( _source ),
				seqindex_( _seqindex ),
				rvs_( _rvs )
		{}

		std::vector< char > const & sequence() const { return sequence_; }
		float score() const { return score_; }
		std::string const & source() const { return source_; }
		unsigned seqindex() const { return seqindex_; }
		bool rvs() const { return rvs_; }

		void print( std::ostream & out = std::cout ) const;

	private:
		std::vector< char > sequence_;
		float score_;
		std::string source_;
		unsigned seqindex_;
		bool rvs_;
};

std::ostream & operator << ( std::ostream & out, Hit const & hit );

////////////////////////////////////////////////////////////////////////////////////////////////////
class HitManager {

	public:
		HitManager() : full_(false), maxhits_(0), outputlevel_(NORMAL) {}

		void full( bool value ) { full_ = value; }
		void maxhits( unsigned value ) { maxhits_ = value; }
		void outputlevel( OutputLevel level ) { outputlevel_ = level; }

		std::list< Hit > const & hits() const { return hits_; }
		bool full() const { return full_; }
		unsigned maxhits() const { return maxhits_; }

		void
		add_hit(
			float score,
			std::vector< char > const & hitseq,
			std::string const & name,
			unsigned seqindex,
			bool rvs = false
		);

		float worst() const { return hits_.front().score(); }
		void print( std::ostream & out = std::cout ) const;

	private:
		std::list< Hit > hits_;
		bool full_;
		unsigned maxhits_;
		OutputLevel outputlevel_;
};

#endif
