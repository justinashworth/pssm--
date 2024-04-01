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
// re-use local copies of these list results for any performance-intensive purposes
std::list<char> nucleotides(){
	std::list<char> nucleotides;
	nucleotides.push_back('A');
	nucleotides.push_back('C');
	nucleotides.push_back('G');
	nucleotides.push_back('T');
	return(nucleotides);
}

std::list<char> base_codes(){
	std::list<char> base_codes;
	base_codes.push_back('A');
	base_codes.push_back('C');
	base_codes.push_back('G');
	base_codes.push_back('T');
	base_codes.push_back('R');
	base_codes.push_back('Y');
	base_codes.push_back('M');
	base_codes.push_back('K');
	base_codes.push_back('S');
	base_codes.push_back('W');
	base_codes.push_back('B');
	base_codes.push_back('D');
	base_codes.push_back('H');
	base_codes.push_back('V');
	base_codes.push_back('N');
	return(base_codes);
}

// don't call this function in any performance-critical operations: use a more efficient long-lived/complex lookup container
std::list<char> degen(char letter){
	std::list<char> degen;
	if(letter=='A'){ degen.push_back('A'); }
	else if(letter=='C'){ degen.push_back('C'); }
	else if(letter=='G'){ degen.push_back('G'); }
	else if(letter=='T'){ degen.push_back('T'); }
	else if(letter=='R'){ degen.push_back('A'); degen.push_back('G'); }
	else if(letter=='Y'){ degen.push_back('C'); degen.push_back('T'); }
	else if(letter=='M'){ degen.push_back('A'); degen.push_back('C'); }
	else if(letter=='K'){ degen.push_back('G'); degen.push_back('T'); }
	else if(letter=='S'){ degen.push_back('C'); degen.push_back('G'); }
	else if(letter=='W'){ degen.push_back('A'); degen.push_back('T'); }
	else if(letter=='B'){
		degen.push_back('C');
		degen.push_back('G');
		degen.push_back('T');
	}
	else if(letter=='D'){
		degen.push_back('A');
		degen.push_back('G');
		degen.push_back('T');
	}
	else if(letter=='H'){
		degen.push_back('A');
		degen.push_back('C');
		degen.push_back('T');
	}
	else if(letter=='V'){
		degen.push_back('A');
		degen.push_back('C');
		degen.push_back('G');
	}
	else if(letter=='N'){
		degen.push_back('A');
		degen.push_back('C');
		degen.push_back('G');
		degen.push_back('T');
	}
	else std::cerr << "Warning: letter " << letter << " not recognized" << std::endl;
	return(degen);
}

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
		case 'N': return true;
		case 'n': return true;
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
		case 'R': return 'R';
		case 'Y': return 'Y';
		case 'M': return 'M';
		case 'K': return 'K';
		case 'S': return 'S';
		case 'W': return 'W';
		case 'B': return 'B';
		case 'D': return 'D';
		case 'H': return 'H';
		case 'V': return 'V';
		case 'N': return 'N';
		case 'n': return 'N';
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
		case 'R': return 'r';
		case 'Y': return 'y';
		case 'M': return 'm';
		case 'K': return 'k';
		case 'S': return 's';
		case 'W': return 'w';
		case 'B': return 'b';
		case 'D': return 'd';
		case 'H': return 'h';
		case 'V': return 'v';
		case 'N': return 'n';
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
  	case 'R': return 'Y';
  	case 'Y': return 'R';
  	case 'S': return 'W';
  	case 'W': return 'S';
  	case 'K': return 'M';
  	case 'M': return 'K';
  	case 'B': return 'V';
  	case 'D': return 'H';
  	case 'H': return 'D';
  	case 'V': return 'B';
		case 'N': return 'N';
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

