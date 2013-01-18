/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/
/** @mainpage
 *  TODO: table of examples with columns:
 *  - adv. algorithm
 *  - Jacobian
 *  - number of dimensions
 *  - equation set
 */

// code licensed under the terms of GNU GPL v3
// copyright holder: University of Warsaw

#pragma once

struct hlf_t {} h;

inline rng_t operator+(
  const rng_t &i, const hlf_t &
) { 
  return i; 
} 

inline rng_t operator-(
  const rng_t &i, const hlf_t &
) { 
  return i-1; 
}

template<class n_t>
inline rng_t operator^(
  const rng_t &r, const n_t &n
) { 
  return rng_t(
    (r - n).first(), 
    (r + n).last()
  ); 
} 