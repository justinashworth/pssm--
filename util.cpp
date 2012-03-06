////////////////////////////////////////////////////////////////////////////////////////////////////
// Justin Ashworth 2007
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdlib> // exit, EXIT_FAILURE

#include "util.h"

std::ostream & operator << (
	std::ostream & ost,
	std::vector< char > const & v
)
{
	for ( std::vector< char >::const_iterator i( v.begin() ), e( v.end() ); i != e; ++i ) ost << *i;
	return ost;
}

////////////////////////////////////////////////////////////////////////////////
bool isnuc( char letter )
{
	switch( letter ) {
		case 'A': return true;
		case 'C': return true;
		case 'G': return true;
		case 'T': return true;
		case 'a': return true;
		case 'c': return true;
		case 'g': return true;
		case 't': return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// probably faster than toupper
char upper( char nucleotide )
{
	switch( nucleotide ) {
		case 'a': return 'A';
		case 'c': return 'C';
		case 'g': return 'G';
		case 't': return 'T';
		case 'A': return 'A';
		case 'C': return 'C';
		case 'G': return 'G';
		case 'T': return 'T';
	}
	std::cerr << "Bad nucleotide " << nucleotide << std::endl;
	exit(EXIT_FAILURE);
	return toupper( nucleotide );
}

////////////////////////////////////////////////////////////////////////////////
// probably faster than tolower
char lower( char nucleotide )
{
	switch( nucleotide ) {
		case 'a': return 'a';
		case 'c': return 'c';
		case 'g': return 'g';
		case 't': return 't';
		case 'A': return 'a';
		case 'C': return 'c';
		case 'G': return 'g';
		case 'T': return 't';
	}
	std::cerr << "Bad nucleotide " << nucleotide << std::endl;
	exit(EXIT_FAILURE);
	return tolower( nucleotide );
}

////////////////////////////////////////////////////////////////////////////////
char comp( char nucleotide )
{
	switch( nucleotide ) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
	}
	std::cerr << "No complement for bad nucleotide " << nucleotide << std::endl;
	exit(EXIT_FAILURE);
	return 'x';
}

////////////////////////////////////////////////////////////////////////////////
// sorting function (should be templated?)
bool secondfloatdesc(
	std::pair< unsigned, float > const & p1,
	std::pair< unsigned, float > const & p2
)
{
//	float f1( p1.second ), f2 ( p2.second );
	return p2.second < p1.second; // descending order
}

