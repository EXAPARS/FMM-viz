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

#include "LB.hpp"

void origin(vec3D * coords, int nbParticles, vector <vec3D> &v)
{
	Particles p;
	p.copyNewCoordinates(coords, nbParticles);
	p.scale();
	v = p.getGlobalCoordsV();	
}

void LBmortonAsync(string file, int nbParticles, vector <vec3D> &v,
	Node<Particles> *& tree, const decompo & nb1ers, const double & dist, double tol,
	const int & first, const int & last, Gaspi_communicator * gComm)
{
	gComm = new Gaspi_communicator(512, nbParticles);
	Particles p(nbParticles, file, gComm);
	p.scale();
	tree = new Node<Particles>(p);
	LB_Base * LBB3 = new LoadBalancer<Particles, MortonAsyncGASPI> (tree, nb1ers, dist, tol, first, last, gComm);
	LBB3->run();	
	v = tree->getContent().getGlobalCoordsV();	
}

void LBmortonSync(string file, int nbParticles, vector <vec3D> &v,
	Node<Particles> *& tree, const decompo & nb1ers, const double & dist, double tol,
	const int & first, const int & last, Gaspi_communicator * gComm)
{
	gComm = new Gaspi_communicator(512, nbParticles);
	Particles p(nbParticles, file, gComm);
	p.scale();
	tree = new Node<Particles>(p);
	LB_Base * LBB2 = new LoadBalancer<Particles, MortonSyncMPI> (tree, nb1ers, dist, tol, first, last, gComm);
	LBB2->run();	
	v = tree->getContent().getGlobalCoordsV();	
}


void LBoctree(vec3D * coords, int targetHeight, int nbParticles, Node<Particles> *& tree, vector <vec3D> &v)
{
	Particles P;
	P.copyNewCoordinates(coords, nbParticles);	
	P.scale();											  	
	tree = new Node<Particles>(P);
	tree -> recDivideOctreeH(targetHeight);
	v = tree->getContent().getGlobalCoordsV();
//	cout << "OCTREE division --> OK" << endl;
}

void LBspectre(string file, vector <vec3D> &v)
{
	Particles p;
	p.loadCoordinatesASCIIWithoutQuantity(file);	
	p.scale();
	v = p.getGlobalCoordsV();
}

void LBspectreHist(string file, vector <vec3D> &v)
{
	double translate = 1436.3972;
	double coeff = 641016.7549;

	Particles p;
	p.loadCoordinatesASCIIWithoutQuantity(file);	
	p.scaleWithParams(translate, coeff);
	v = p.getGlobalCoordsV();
}
