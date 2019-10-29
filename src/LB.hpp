/*
  Copyright 2015 - UVSQ
  Authors list: Nathalie Möller, Eric Petit

  This file is part of the FMM-viz.

  FMM-viz is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later version.

  FMM-viz is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License along with
  the FMM-viz. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LB_HPP
#define LB_HPP

#include <vector>
#include "Node.hpp"
#include "Particles.hpp"
#include "vec3D.hpp"
#include "Decomposition.hpp"

#include "LoadBalancerBase.hpp"
#include "LoadBalancer.hpp"
#include "LBHistExact.hpp"
#include "LBHistApprox.hpp"
#include "LBMortonSyncMPI.hpp"
#include "LBMortonAsyncGASPI.hpp"



using namespace std;


void origin(vec3D * coords, int nbParticles, vector <vec3D> &v);

template <typename T>
void LBB(vec3D * coords, int nbParticles, vector <vec3D> &v,
	Node<Particles> *& tree, const decompo & nb1ers, const double & dist, double tol,
	const int & first, const int & last, Gaspi_communicator * gComm)
{

	Particles p;
	p.copyNewCoordinates(coords, nbParticles);
	p.scale();												  
	tree = new Node<Particles>(p);
	LB_Base * LBB0 = new LoadBalancer<Particles, T> (tree, nb1ers, dist, tol, first, last, gComm);
	LBB0 -> run();	
	v = tree->getContent().getGlobalCoordsV();
}

void LBmortonAsync(string file, int nbParticles, vector <vec3D> &v,
	Node<Particles> *& tree, const decompo & nb1ers, const double & dist, double tol,
	const int & first, const int & last, Gaspi_communicator * gComm);

void LBoctree(vec3D * coords, int targetHeight, int nbParticles, Node<Particles> *& tree, vector <vec3D> &v);
void LBspectre(string file, vector <vec3D> &v);
void LBspectreHist(string file, vector <vec3D> &v);

	
#endif
