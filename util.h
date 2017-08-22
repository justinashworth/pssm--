////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_util
#define INCLUDED_util

#include <iosfwd>
#include <vector>
#include <list>

enum OutputLevel {
	MINIMAL,
	NORMAL,
	VERBOSE
};

std::ostream & operator << (
	std::ostream & ost,
	std::vector< char > const & v
);

std::list<char> nucleotides();
std::list<char> base_codes();
std::list<char> degen(char);
bool isnuc(char);
char upper(char);
char lower(char);
char comp(char);

bool secondfloatdesc(
	std::pair< unsigned, float > const & p1,
	std::pair< unsigned, float > const & p2
);

#endif
