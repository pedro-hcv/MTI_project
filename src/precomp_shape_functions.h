/*
 * ================================================
 * 					COPYRIGHT:
 * Institute of Machine Tools & Manufacturing (IWF)
 * Department of Mechanical & Process Engineering
 * 					ETH ZURICH
 * ================================================
 *
 *  This file is part of "mfree_iwf-ul-cut-refine".
 *
 * 	mfree_iwf is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	mfree_iwf-ul-cut-refine is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *  along with mfree_iwf-ul-cut-refine.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *  This is the source code used to produce the results
 *  of the metal cutting simulation presented in:
 *
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  "Meshfree Simulation of Metal Cutting:
 *  An Updated Lagrangian Approach with Dynamic Refinement"
 *
 * 	Authored by:
 * 	Mohamadreza Afrasiabi
 * 	Dr. Matthias Roethlin
 * 	Hagen Klippel
 * 	Prof. Dr. Konrad Wegener
 *
 * 	Published by:
 * 	International Journal of Mechanical Sciences
 * 	28 June 2019
 * 	https://doi.org/10.1016/j.ijmecsci.2019.06.045
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * 	For further descriptions, you may refer to the manuscript
 * 	or the previous works of the same research group
 * 	at IWF, ETH Zurich.
 *
 */

#ifndef PRECOMPSHAPEFUNCTIONS_H_
#define PRECOMPSHAPEFUNCTIONS_H_

#include "particle.h"
#include "kernel.h"

#include <assert.h>
#include <vector>
#include <glm/glm.hpp>

/*
  This implementation intends to compute the smoothing kernel (W)
  and its 1st spatial derivatives using CSPM.

  In summary:
  ==========================
  1- precomp_sph: implementation of Eq. (22) from the paper.
  	  	  	  	  calculate and assign the kernel values of
  	  	  	  	  each particle "i" with respect to its neighbors "j" with standard SPH


  2- precomp_cspm: implementation of Eqs. (28)-(29) form the paper.
  	  	  	  	   calculate and assign the corrected kernel and its 1st derivatives values of
  	  	  	  	   each particle "i" with respect to its neighbors "j" with standard CSPM
*/

void precomp_sph(std::vector<particle> &particles,  unsigned int n);
void precomp_cspm(std::vector<particle> &particles, unsigned int n);

#endif /* PRECOMPSHAPEFUNCTIONS_H_ */
