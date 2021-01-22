/*********************************************************************************

			   Main View for Point Cloud Remeshing and Mesh 
			                Reconstruction

						Updating in 2020/09/09

						   By Dr. Chenlei Lv

			The functions includes:
			1. The pipeline of point cloud-based shape reconstruction;
			


*********************************************************************************/
#pragma once
#include "View.h"
#include <stdio.h>
#include <stdlib.h>
#include "trackball.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Method_AIVS_RNPro.hpp"
using namespace std;

int main()
{
	cout << "start!" << endl;

	//input parameters
	//support form: obj, xyz, ply, txt, off
	string pathObj = "Bunny";	
	string pathFile = "data//stand//" + pathObj + ".ply"; 
	int pointNum = 10000;//simplification Number
	int classificationIndex = 1;//1: isotropic; 2: Curvature sensitive; 3: Edge keeping. 
	double multi = 10;//default multi = 2 ~ 10， used for error mesh remove.
	int optIter = 5;//default optIter = 5, geometric optimize iteration number.
	
	//resampling and remeshing
	AIVS_RNPro ar;
	ar.AIVS_RNPro_init(pathFile);		
	ar.AIVS_RNPro_Remesh(pathObj, pointNum, classificationIndex, multi, optIter);
	
  	cout << "Finish!" << endl;

}


















