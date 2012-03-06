////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iomanip> // std::showpoint, std::fixed, std::setprecision
#include <algorithm> // std::reverse, etc.

#include "Hits.h"
#include "util.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
void
Hit::print( std::ostream & out ) const
{
	out << std::showpoint << std::fixed << std::setprecision(2) << score_
	    << " " << sequence_ << " " << source_ << " " << seqindex_;
	if ( rvs_ ) out << " (rvs)";
	out << std::endl;
}

std::ostream & operator << ( std::ostream & out, Hit const & hit )
{
	hit.print( out );
	return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
HitManager::add_hit(
	float score,
	std::vector< char > const & hitseq,
	std::string const & name,
	unsigned seqindex,
	bool rvs
)
{

	std::list< Hit >::iterator insert_itr( hits_.end() );
	// first hit, or better than current best
	if ( hits_.size() == 0 || score < hits_.back().score() ) insert_itr = hits_.end();
	// in between: insert in proper order
	else {
	// this could be made MUCH faster for large lists, as a binary search that assumes strict ordering
	// dilemma: list is constant insert time, but no direct access to elements (must traverse list)
		for ( std::list< Hit >::iterator hit_itr( hits_.begin() );
		      hit_itr != hits_.end(); ++hit_itr ) {
			// found the first superior hit_itr, insert before it
			if ( score >= hit_itr->score() ) { insert_itr = hit_itr; break; }
		}
	}

	if ( !rvs ) {
		hits_.insert( insert_itr, Hit( hitseq, score, name, seqindex ) );
	} else {
		std::vector<char> rvsseq( hitseq );
		std::reverse( rvsseq.begin(), rvsseq.end() );
		std::transform( rvsseq.begin(), rvsseq.end(), rvsseq.begin(), comp );
		hits_.insert( insert_itr, Hit( rvsseq, score, name, seqindex, true ) );
	}

	if ( hits_.size() > maxhits_ ) {
		full_ = true;
		hits_.pop_front();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
HitManager::print( std::ostream & out ) const
{
	for ( std::list< Hit >::const_iterator h( hits_.begin() ), e( hits_.end() ); h != e; ++h ) out << *h;
}

