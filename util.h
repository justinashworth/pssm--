////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_util
#define INCLUDED_util

#include <iosfwd>
#include <vector>

enum OutputLevel {
	MINIMAL,
	NORMAL,
	VERBOSE
};

std::ostream & operator << (
	std::ostream & ost,
	std::vector< char > const & v
);

bool isnuc( char letter );
char upper( char nucleotide );
char lower( char nucleotide );
char comp( char nucleotide );

bool secondfloatdesc(
	std::pair< unsigned, float > const & p1,
	std::pair< unsigned, float > const & p2
);

#endif
