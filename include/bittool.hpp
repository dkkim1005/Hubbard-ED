#pragma oncce

// last updated : 19 APR 2012
//      * parity_count added 
//      * bitcount updated

/*
// table of the number of bits for the numbers 0..255
static const unsigned numbits_lookup_table[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2,
    3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
    3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,
    4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
    3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5,
    6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
    4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
    6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3,
    4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
    6, 7, 6, 7, 7, 8
}; 

static const unsigned nCm[21][21] = { {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,5,10,10,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,6,15,20,15,6,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,7,21,35,35,21,7,1,0,0,0,0,0,0,0,0,0,0,0,0,0}, {1,8,28,56,70,56,28,8,1,0,0,0,0,0,0,0,0,0,0,0,0}, {1,9,36,84,126,126,84,36,9,1,0,0,0,0,0,0,0,0,0,0,0}, {1,10,45,120,210,252,210,120,45,10,1,0,0,0,0,0,0,0,0,0,0}, {1,11,55,165,330,462,462,330,165,55,11,1,0,0,0,0,0,0,0,0,0}, {1,12,66,220,495,792,924,792,495,220,66,12,1,0,0,0,0,0,0,0,0}, {1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1,0,0,0,0,0,0,0}, {1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1,0,0,0,0,0,0}, {1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365,455,105,15,1,0,0,0,0,0}, {1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1,0,0,0,0}, {1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,136,17,1,0,0,0}, {1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824,18564,8568,3060,816,153,18,1,0,0}, {1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1,0}, {1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1} };

// calculate the number of 1-bits : works for both of 32/64-bit machines
inline unsigned bitcount(const unsigned& s) {
	unsigned count = 0;
	unsigned n=s;
	while (n) {
		count += numbits_lookup_table[n & 0xffu];
		n >>= 8;
	}
	return count;
}

*/

// FROM HERE
// REF: see "Bit Twiddling Hacks" by Sean Eron Anderson
// http://graphics.stanford.edu/~seander/bithacks.html

namespace bittool
{
static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
	    B6(0), B6(1), B6(1), B6(2)
};

static const bool ParityTable256[256] = 
{
#   define P2(n) n, n^1, n^1, n
#   define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#   define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
	    P6(0), P6(1), P6(1), P6(0)
};

inline unsigned bitcount(const unsigned & s)
{
  return BitsSetTable256[s & 0xff] + 
    BitsSetTable256[(s >> 8) & 0xff] + 
    BitsSetTable256[(s >> 16) & 0xff] + 
    BitsSetTable256[s >> 24];
}

inline bool pairty_count_all(const signed & s)
{
  unsigned v = s;
  v ^= v >> 16;
  v ^= v >> 8;
  return ParityTable256[v & 0xff];
}

inline bool parity_count(const unsigned & i, const unsigned & s)
{
  unsigned v = s & ((0x1u << i) - 0x1u);
  v ^= v >> 16;
  v ^= v >> 8;
  return ParityTable256[v & 0xff];
}

// TO HERE

// Retrieve the bit (0 or 1) at i-th position
inline unsigned bitget(const unsigned & i, const unsigned & s) { return (s & 0x1u<<i)>>i; }

// Flip the bit at i-th position
inline unsigned bitflip(const unsigned & i, const unsigned & s) { return s ^ 0x1u<<i; }

// Hacker's Delight (Addison-Wesley,2003), Ch. 2, see hackersdelight.org
// Next higher number having the same total number of 1-bits
inline unsigned next_set_of_n_elements(const unsigned & x)
{
  unsigned smallest, ripple, ones;
  if (x == 0)
    return 0;
  smallest = x & -x;
  ripple = x + smallest;
  ones = x ^ ripple;
  ones = (ones >> 2)/smallest;
  return ripple | ones;
}
} // end namespace bittool
