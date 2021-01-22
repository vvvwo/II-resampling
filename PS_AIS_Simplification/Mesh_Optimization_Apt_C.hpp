/*********************************************************************************

	  Mesh Optimization to achieve adaptive isotropic remeshing (Curvature)

						 Updating in 2020/12/01

						    By Dr. Chenlei Lv

			The functions includes:
			[Dunyach 2013] Adaptive Remeshing for Real-Time Mesh Deformation.

*********************************************************************************/

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h> 
#include <algorithm>
#include <pcl/io/pcd_io.h>  //File input$output
#include <pcl/octree/octree_search.h>  //octree define
#include <pcl/point_types.h>  //point type
#include <pcl/kdtree/kdtree_flann.h>
#include "Mesh_Geometric.hpp"
#define PI 3.1415926
using namespace std;

class MeshOptimization_Apt_C {

private:

	int originalNumber;
	string fileName;
	vector<vector<double>> points;//point cloud
	//vector<double> pL;
	vector<double> cu_ave;
	vector<vector<int>> faceInfor;//trangulars
	vector<vector<int>> pointNeighbor;//neighbor faces
	int optIter;
	//double L_ave; //the ave edge length of all borders
	double L_regular;	
	vector<double> L_range;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeG;
	int iterStep = 10;

public:

	void MeshOptimization_Apt_C_init(string fileN, vector<vector<double>> pointsi, vector<vector<int>> faceInfori) {
		originalNumber = pointsi.size();
		fileName = fileN;
		cout << "Mesh optimization start!" << endl;
		points = pointsi;
		faceInfor = faceInfori;//the index of first point is 1		
		//init neighbor structure, correct point index to first point index is 0
		pointNeighbor.resize(points.size());
		for (int i = 0; i < faceInfor.size(); i++) {
			//faceInfor[i][0]--;
			//faceInfor[i][1]--;
			//faceInfor[i][2]--;
			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];
			pointNeighbor[b1].push_back(b2);
			pointNeighbor[b1].push_back(b3);
			pointNeighbor[b2].push_back(b3);
			pointNeighbor[b2].push_back(b1);
			pointNeighbor[b3].push_back(b1);
			pointNeighbor[b3].push_back(b2);
		}
		//init ave length of trangular border
		MeshOptimization_Apt_AveEdge_init();
		
		MeshOptimization_Apt_L_init();

		//Kd-tree construct
		kdtreeG = MeshOptimization_Apt_KdTreeConstruct(pointsi);

		cout << "L_regular:" << L_regular << endl;
	}

	void MeshOptimization_Apt_C_Start(int iter, bool judge) {

		optIter = iter;

		cout << "Mesh Split_Collapse start!" << endl;
		while (iter) {
			cout << "iter:" << iter << endl;
			//processing split, collapse, Flip and Tangent smoothing
			cout << "Split start:";
			MeshOptimization_Apt_Split();

			cout << "Collapse start:";
			MeshOptimization_Apt_Collapse(judge);

			cout << "Flip start:";
			MeshOptimization_Apt_Flip();

			cout << "Tangent smoothing start:";
			MeshOptimization_Apt_TangentSmoothing();

			iter--;
		}
		//cout << "tangent moving:" << endl;
		//MeshOptimization_TangentSmoothing();
		MeshOptimization_Apt_UpdateStructure_Face();
		//cout << "Mesh Save." << endl;
		//MeshOptimization_SaveOBJ_Face();
		cout << endl;

	}

	vector<vector<double>> MeshOptimization_Apt_C_GetPoints() {
		return points;
	}

	vector<vector<int>> MeshOptimization_Apt_C_GetFaces() {
		return faceInfor;
	}

private:

	//construct kd-tree and return adaptive L
	pcl::KdTreeFLANN<pcl::PointXYZ>  MeshOptimization_Apt_KdTreeConstruct(vector<vector<double>> seedpoints) {

		pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeSeed;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->width = seedpoints.size();
		cloud->height = 1;
		cloud->points.resize(cloud->width * cloud->height);
		// fills a PointCloud with random data
		for (int i = 0; i < seedpoints.size(); i++)
		{
			pcl::PointXYZ pxyz;
			cloud->points[i].x = seedpoints[i][0];
			cloud->points[i].y = seedpoints[i][1];
			cloud->points[i].z = seedpoints[i][2];

		}
		kdtreeSeed.setInputCloud(cloud);
		return kdtreeSeed;
	}

	double MeshOptimization_Apt_KdTreeSearch(vector<double> pointP) {

		vector<int> pointIdxNKNSearch;
		vector<float> pointNKNSquaredDistance;
		int kn = 1;
		pcl::PointXYZ searchPoint;
		searchPoint.x = pointP[0];
		searchPoint.y = pointP[1];
		searchPoint.z = pointP[2];
		kdtreeG.nearestKSearch(searchPoint, kn, pointIdxNKNSearch, pointNKNSquaredDistance);
		int L_index = pointIdxNKNSearch[0];
		return cu_ave[L_index];

	}

	double MeshOptimization_Apt_KdTreeSearch(int pindex) {

		vector<int> pointIdxNKNSearch;
		vector<float> pointNKNSquaredDistance;
		int kn = 1;
		pcl::PointXYZ searchPoint;
		searchPoint.x = points[pindex][0];
		searchPoint.y = points[pindex][1];
		searchPoint.z = points[pindex][2];
		kdtreeG.nearestKSearch(searchPoint, kn, pointIdxNKNSearch, pointNKNSquaredDistance);
		int L_index = pointIdxNKNSearch[0];
		return cu_ave[L_index];

	}

	//Split, return the sum of point number
	void MeshOptimization_Apt_Split() {

		int iter = 0;
		while (1) {

			int Split_list_index = 0;
			int pointsEndIndex = points.size();//the index for last point index, by update point list
			int numberProcess = 0;

			vector<vector<int>> Split_list;
			for (int i = 0; i < pointNeighbor.size(); i++) {
				vector<int> b2_list;
				int b1 = i;
				for (int j = 0; j < pointNeighbor[i].size(); j++) {
					int b2 = pointNeighbor[i][j];
					if (!MeshOptimization_Apt_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_Apt_EdgeLength(b1, b2);
						double L_j1 = MeshOptimization_Apt_KdTreeSearch(b1);
						double L_j2 = MeshOptimization_Apt_KdTreeSearch(b2);
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (length12 >= 1.33 * L_j1) {
							vector<int> edgeij;
							edgeij.push_back(b1);
							edgeij.push_back(b2);
							Split_list.push_back(edgeij);
						}
					}
				}
			}

			cout << Split_list.size();
			if (Split_list.size() == 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}

			//Split processing
			for (int i = 0; i < Split_list.size(); i++) {
				int bs1 = Split_list[i][0];
				int bs2 = Split_list[i][1];
				//Split: 1 achieve middle point
				vector<double> middlePoint(3);
				middlePoint[0] = (points[bs1][0] + points[bs2][0]) / 2;
				middlePoint[1] = (points[bs1][1] + points[bs2][1]) / 2;
				middlePoint[2] = (points[bs1][2] + points[bs2][2]) / 2;
				points.push_back(middlePoint);

				vector<int> pointNeighbor_b1 = pointNeighbor[bs1];
				vector<int> pointNeighbor_middle;//store the structure of middle point

				//Split: 2 store related points
				vector<int> bs3_list;
				bs3_list.push_back(bs1);
				bs3_list.push_back(bs2);
				for (int i = 0; i < pointNeighbor_b1.size() / 2; i++) {
					int b1i1 = pointNeighbor_b1[2 * i];
					int b1i2 = pointNeighbor_b1[2 * i + 1];
					if (b1i1 == bs2) {
						bs3_list.push_back(b1i2);
					}
					else if (b1i2 == bs2) {
						bs3_list.push_back(b1i1);
					}
					else {
						continue;
					}
				}

				//Split: 3 update the structure			
				vector<vector<int>> T_n;
				for (int i = 0; i < bs3_list.size(); i++) {
					int b_index = bs3_list[i];
					if (b_index == bs1 || b_index == bs2) {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if (b_index_j1 == bs1 || b_index_j1 == bs2) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
							}
							if (b_index_j2 == bs1 || b_index_j2 == bs2) {
								pointNeighbor[b_index][2 * j + 1] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
					else {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if ((b_index_j1 == bs1 || b_index_j1 == bs2) && (b_index_j2 == bs1 || b_index_j2 == bs2)) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor[b_index].push_back(b_index_j1);
								pointNeighbor[b_index].push_back(pointsEndIndex);
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
				}
				//Split: 4 update neighbor structure
				pointNeighbor.push_back(MeshOptimization_Apt_Structure_RemoveRepeat(pointNeighbor_middle));
				//pointNeighbor.push_back(pointNeighbor_middle);
				pointsEndIndex++;
			}

			//Remove Repeat from neighbor
			for (int i = originalNumber; i < pointNeighbor.size(); i++) {
				vector<int> v = pointNeighbor[i];
				vector<int> v_RemoveRepeat = MeshOptimization_Apt_Structure_RemoveRepeat(v);
				pointNeighbor[i].clear();
				pointNeighbor[i] = v_RemoveRepeat;
			}
			iter++;
			if (iter > iterStep) {
				break;			
			}
		}
		cout << endl;
	}

	//Collapse, return the point number which has been deleted from the original point cloud
	void MeshOptimization_Apt_Collapse(bool judge) {

		int CollapseNum = points.size() - originalNumber;

		while (1) {

			if (CollapseNum == 0 && judge) {
				break;//control the certain point number			
			}

			int Collapse_Number = 0;

			vector<bool> pointActive(points.size(), true);//state of point, if ture, can be Split-Collapse, alse can not		
			vector<int> point_Delete;//store the index of delete points.
			int pointsEndIndex = points.size();

			//store all edges satisfy the condition of Split-Collapse		
			vector<vector<int>> Collapse_List;

			//cout << "Achieve all edges for Collapse" << endl;
			//achieve all edges for Collapse
			for (int i = 0; i < pointNeighbor.size(); i++) {

				vector<int> b2_list;
				int b1 = i;
				for (int j = 0; j < pointNeighbor[i].size(); j++) {
					int b2 = pointNeighbor[i][j];
					if (!MeshOptimization_Apt_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_Apt_EdgeLength(b1, b2);
						double L_j1 = MeshOptimization_Apt_KdTreeSearch(b1);
						double L_j2 = MeshOptimization_Apt_KdTreeSearch(b2);
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (L_j1 > L_j2) {
							L_j1 = L_j2;
						}
						if (length12 < 0.8 * L_j1) {
							//Test if there have invaild condition
							int valence = MeshOptimization_Apt_Edge_Valence(b1, b2);
							if (valence <= 2) {
								vector<int> edgeij;
								edgeij.push_back(b1);
								edgeij.push_back(b2);
								Collapse_List.push_back(edgeij);
							}
						}
					}
				}
			}

			for (int i = 0; i < Collapse_List.size(); i++) {
				//*****************************Collapse test**********************************			
				int bc1 = -1;
				int bc2 = -1;
				int bc1_i = Collapse_List[i][0];
				int bc2_i = Collapse_List[i][1];
				int valenceE = MeshOptimization_Apt_Edge_Valence(bc1_i, bc2_i);
				if (pointActive[bc1_i] && pointActive[bc2_i] && valenceE <= 2) {
					bc1 = bc1_i;
					bc2 = bc2_i;
				}
				else {
					continue;
				}

				vector<int> pointNeighbor_bc1 = pointNeighbor[bc1];
				vector<int> pointNeighbor_middlec;//store the structure of new point
				vector<double> middlePointc(3);
				middlePointc[0] = (points[bc1][0] + points[bc2][0]) / 2;
				middlePointc[1] = (points[bc1][1] + points[bc2][1]) / 2;
				middlePointc[2] = (points[bc1][2] + points[bc2][2]) / 2;

				//Collapse test
				//add the related points in the list
				vector<int> pc_list;
				//bc1 related point
				for (int i = 0; i < pointNeighbor[bc1].size(); i++) {
					int p_i = pointNeighbor[bc1][i];
					if (p_i == bc1 || p_i == bc2) {
						continue;
					}
					if (!MeshOptimization_Apt_Exist(p_i, pc_list) && pointActive[p_i]) {
						pc_list.push_back(p_i);
					}
				}
				//bc2 related point
				for (int i = 0; i < pointNeighbor[bc2].size(); i++) {
					int p_i = pointNeighbor[bc2][i];
					if (p_i == bc1 || p_i == bc2) {
						continue;
					}
					if (!MeshOptimization_Apt_Exist(p_i, pc_list) && pointActive[p_i]) {
						pc_list.push_back(p_i);
					}
				}

				//check the long edge border length should not longer than 4/3
				bool Collapse_Invaild = false;
				for (int i = 0; i < pc_list.size(); i++) {
					vector<double> pN2 = points[pc_list[i]];
					double radius12 = sqrt((middlePointc[0] - pN2[0]) * (middlePointc[0] - pN2[0]) +
						(middlePointc[1] - pN2[1]) * (middlePointc[1] - pN2[1]) +
						(middlePointc[2] - pN2[2]) * (middlePointc[2] - pN2[2]));

					double L_j1 = MeshOptimization_Apt_KdTreeSearch(middlePointc);
					double L_j2 = MeshOptimization_Apt_KdTreeSearch(pN2);
					if (L_j1 > L_j2) {
						L_j1 = L_j2;
					}

					if (radius12 > 1.33 * L_j1) {
						Collapse_Invaild = true;
						break;
					}
				}
				if (Collapse_Invaild) {
					continue;
				}

				Collapse_Number++;
				//Collapse processing;
				//Collapse: 1 init
				point_Delete.push_back(bc1);
				point_Delete.push_back(bc2);
				pointActive[bc1] = false;
				pointActive[bc2] = false;
				points.push_back(middlePointc);
				pointActive.push_back(true);

				//Collapse: 2 update the structure
				for (int i = 0; i < pc_list.size(); i++) {
					int p_index = pc_list[i];
					bool rejudge = false;

					for (int j = 0; j < pointNeighbor[p_index].size() / 2; j++) {
						int pc1 = pointNeighbor[p_index][2 * j];
						int pc2 = pointNeighbor[p_index][2 * j + 1];

						if ((pc1 == bc1 || pc1 == bc2) && (pc2 == bc1 || pc2 == bc2)) {
							pointNeighbor[p_index][2 * j] = -1;
							pointNeighbor[p_index][2 * j + 1] = -1;
							rejudge = true;
						}
						else if (pc1 == bc1 || pc1 == bc2) {
							pointNeighbor[p_index][2 * j] = pointsEndIndex;
							pointNeighbor_middlec.push_back(pc2);
							pointNeighbor_middlec.push_back(p_index);
						}
						else if (pc2 == bc1 || pc2 == bc2) {
							pointNeighbor[p_index][2 * j + 1] = pointsEndIndex;
							pointNeighbor_middlec.push_back(p_index);
							pointNeighbor_middlec.push_back(pc1);
						}
						else {
							continue;
						}
					}
					if (rejudge) {
						vector<int> pointNeighbor_new;
						for (int j = 0; j < pointNeighbor[p_index].size(); j++) {
							if (pointNeighbor[p_index][j] != -1) {
								pointNeighbor_new.push_back(pointNeighbor[p_index][j]);
							}
						}
						pointNeighbor[p_index].clear();
						pointNeighbor[p_index] = pointNeighbor_new;
					}
				}
				pointNeighbor[bc1].clear();
				pointNeighbor[bc1].push_back(-1);
				pointNeighbor[bc1].push_back(-1);

				pointNeighbor[bc2].clear();
				pointNeighbor[bc2].push_back(-1);
				pointNeighbor[bc2].push_back(-1);

				pointNeighbor.push_back(MeshOptimization_Apt_Structure_RemoveRepeat(pointNeighbor_middlec));

				pointsEndIndex++;
				CollapseNum--;
				if (CollapseNum == 0) {
					break;//control the certain point number			
				}
			}

			cout << Collapse_Number;
			if (Collapse_Number <= 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}
			//Remove Repeat from neighbor
			for (int i = originalNumber; i < pointNeighbor.size(); i++) {
				vector<int> v = pointNeighbor[i];
				vector<int> v_RemoveRepeat = MeshOptimization_Apt_Structure_RemoveRepeat(v);
				pointNeighbor[i].clear();
				pointNeighbor[i] = v_RemoveRepeat;
			}
			MeshOptimization_Apt_UpdateStructure(point_Delete);
		}
		cout << endl;
	}

	//Flip
	void MeshOptimization_Apt_Flip() {
		while (1) {
			int FilpNum = 0;
			for (int i = 0; i < points.size(); i++) {
				//Flip test					
				int b1 = i;
				//cout << b1;			
				//Flip start
				vector<int> pNb1;//add points into list
				for (int j = 0; j < pointNeighbor[b1].size(); j++) {
					int b2_j = pointNeighbor[b1][j];
					if (!MeshOptimization_Apt_Exist(b2_j, pNb1)) {
						pNb1.push_back(b2_j);
					}
				}

				for (int j = 0; j < pNb1.size(); j++) {
					int b2 = pNb1[j];
					if (b2 < b1) {
						continue;
					}
					//cout << b2;
					vector<int> b3List;//store the other two points for flip
					for (int k = 0; k < pointNeighbor[b1].size() / 2; k++) {
						int b11 = pointNeighbor[b1][2 * k];
						int b12 = pointNeighbor[b1][2 * k + 1];
						if (b11 == b2 && !MeshOptimization_Apt_Exist(b12, b3List)) {
							b3List.push_back(b12);
						}
						if (b12 == b2 && !MeshOptimization_Apt_Exist(b11, b3List)) {
							b3List.push_back(b11);
						}
					}
					if (b3List.size() == 2) {
						//Flip processing
						//compute 4 point valence
						int b3 = b3List[0];
						int b4 = b3List[1];
						//cout <<"Flip:"<< b1 << "," << b2 << "," << b3 << "," << b4 << endl;
						if (MeshOptimization_Apt_If_Flip(b1, b2, b3, b4)) {
							//Flip update
							//cout << b1 << "," << b2 << "," << b3 << "," << b4 << endl;
							MeshOptimization_Apt_Flip_Processing(b1, b2, b3, b4);
							FilpNum++;
						}
					}
					else {
						continue;
					}
				}
			}
			cout << FilpNum;
			if (FilpNum == 0) {
				cout << ".";
				break;
			}
			else {
				cout << ",";
			}
		}
		cout << endl;
	}

	//Update structure
	void MeshOptimization_Apt_UpdateStructure(vector<int> pointDelete) {

		vector<int> pointsIndex(points.size());
		vector<bool> pointjudge(points.size(), true);
		for (int i = 0; i < pointsIndex.size(); i++) {
			pointsIndex[i] = i;
		}

		for (int i = 0; i < pointDelete.size(); i++) {
			int p_index = pointDelete[i];
			pointjudge[p_index] = false;
			pointsIndex[p_index] = -1;
			for (int j = p_index + 1; j < pointsIndex.size(); j++) {
				pointsIndex[j]--;
			}
		}

		//update point list
		vector<vector<double>> pointNew;
		vector<vector<int>> pointNeighborNew;
		vector<vector<int>> faceInforNew;
		for (int i = 0; i < pointjudge.size(); i++) {
			if (pointjudge[i]) {
				pointNew.push_back(points[i]);
				vector<int> pointNeighborNew_i = pointNeighbor[i];
				for (int j = 0; j < pointNeighborNew_i.size(); j++) {
					pointNeighborNew_i[j] = pointsIndex[pointNeighborNew_i[j]];
					//int index = pointNeighbor[i][j];
					//pointNeighborNew_i.push_back(pointsIndex[index]);
				}
				pointNeighborNew.push_back(pointNeighborNew_i);
			}
		}
		points.clear();
		points = pointNew;
		pointNeighbor.clear();
		pointNeighbor = pointNeighborNew;
	}

	//Finally output the 
	void MeshOptimization_Apt_UpdateStructure_Face() {

		vector<vector<int>> pointNeighborNew = pointNeighbor;
		vector<vector<int>> faceInforNew;
		//update face list
		for (int i = 0; i < pointNeighborNew.size(); i++) {
			int b1 = i;
			vector<int> b1_N = pointNeighborNew[b1];
			for (int j = 0; j < pointNeighborNew[i].size() / 2; j++) {
				int b2 = pointNeighborNew[i][2 * j];
				int b3 = pointNeighborNew[i][2 * j + 1];
				if (b2 == -1 || b3 == -1) {
					continue;
				}
				else {
					//f1 << "f " << i + 1 << " " << pointNeighborNew[i][2 * j] + 1 << " " << pointNeighborNew[i][2 * j + 1] + 1 << endl;
					vector<int> faceInforNew_i;
					faceInforNew_i.push_back(i);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j]);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j + 1]);
					faceInforNew.push_back(faceInforNew_i);
					for (int k = 0; k < pointNeighborNew[b2].size() / 2; k++) {
						int b2n2 = pointNeighborNew[b2][2 * k];
						int b2n3 = pointNeighborNew[b2][2 * k + 1];
						if ((b2n2 == b1 || b2n2 == b3) && (b2n3 == b1 || b2n3 == b3)) {
							pointNeighborNew[b2][2 * k] = -1;
							pointNeighborNew[b2][2 * k + 1] = -1;
						}
					}
					for (int k = 0; k < pointNeighborNew[b3].size() / 2; k++) {
						int b3n2 = pointNeighborNew[b3][2 * k];
						int b3n3 = pointNeighborNew[b3][2 * k + 1];
						if ((b3n2 == b1 || b3n2 == b2) && (b3n3 == b1 || b3n3 == b2)) {
							pointNeighborNew[b3][2 * k] = -1;
							pointNeighborNew[b3][2 * k + 1] = -1;
						}
					}
				}
			}
		}
		faceInfor.clear();
		faceInfor = faceInforNew;
	}

	//*******************************************Basic Function*****************************************************

	//compute ave edge length of a mesh	
	void MeshOptimization_Apt_AveEdge_init() {//computer the ave edge l for trangular edge estimation

		//select 1000 trangulars for estimation

		double l_sum = 0;
		for (int i = 0; i < faceInfor.size(); i++) {
			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];
			double l1 = MeshOptimization_Apt_EdgeLength(b1, b2);
			double l2 = MeshOptimization_Apt_EdgeLength(b2, b3);
			double l3 = MeshOptimization_Apt_EdgeLength(b3, b1);
			double l_ave_i = (l1 + l2 + l3) / 3;
			l_sum = l_sum + l_ave_i;
		}
		L_regular = l_sum / faceInfor.size();
		
	}

	//return the length of an edge (b1, b2)
	double MeshOptimization_Apt_EdgeLength(int b1, int b2) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];

		double lengthE = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
			(p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));

		return lengthE;

	}

	//judge the point is in a list
	bool MeshOptimization_Apt_Exist(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return true;
			}
		}
		return false;
	}

	int MeshOptimization_Apt_Exist_Index(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return i;
			}
		}
		return -1;
	}

	//return the Valence value of a edge (p1, p2). 1: edge; 2: normal; 3: invalid
	int MeshOptimization_Apt_Edge_Valence(int p1, int p2) {
		int valence = 0;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			if (pointNeighbor[p1][2 * i] == p2 || pointNeighbor[p1][2 * i + 1] == p2) {
				valence++;
			}
		}
		return valence;
	}

	//return the Valence value of a point
	int MeshOptimization_Apt_Point_Valence(int p1) {


		vector<int> b_num;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			int p2 = pointNeighbor[p1][2 * i];
			int p3 = pointNeighbor[p1][2 * i + 1];
			if (!MeshOptimization_Apt_Exist(p2, b_num)) {
				b_num.push_back(p2);
			}
			if (!MeshOptimization_Apt_Exist(p3, b_num)) {
				b_num.push_back(p3);
			}
		}
		return b_num.size();

	}

	//return if a point is an edge point
	bool MeshOptimization_Apt_Edge_Point(int p1) {

		for (int i = 0; i < pointNeighbor[p1].size(); i++) {
			int p2 = pointNeighbor[p1][i];
			int valenceE = MeshOptimization_Apt_Edge_Valence(p1, p2);
			if (valenceE == 1) {
				return true;
			}
		}
		return false;
	}

	//Flip judge and processing
	bool MeshOptimization_Apt_If_Flip(int b1, int b2, int b3, int b4) {

		double disb34 = MeshOptimization_Apt_EdgeLength(b3, b4);

		double L_j1 = MeshOptimization_Apt_KdTreeSearch(b3);
		double L_j2 = MeshOptimization_Apt_KdTreeSearch(b4);
		if (L_j1 > L_j2) {
			L_j1 = L_j2;
		}

		if (disb34 > 1.33 * L_j1 || disb34 < 0.8 * L_j1) {
			return false;
		}
		bool b1e = MeshOptimization_Apt_Edge_Point(b1);
		bool b2e = MeshOptimization_Apt_Edge_Point(b2);
		bool b3e = MeshOptimization_Apt_Edge_Point(b3);
		bool b4e = MeshOptimization_Apt_Edge_Point(b4);
		int b1ValenceTarget = 6;
		if (b1e) {
			b1ValenceTarget = 4;
		}
		int b2ValenceTarget = 6;
		if (b2e) {
			b2ValenceTarget = 4;
		}
		int b3ValenceTarget = 6;
		if (b3e) {
			b3ValenceTarget = 4;
		}
		int b4ValenceTarget = 6;
		if (b4e) {
			b4ValenceTarget = 4;
		}
		int b1ValenceBefore = MeshOptimization_Apt_Point_Valence(b1);
		int b2ValenceBefore = MeshOptimization_Apt_Point_Valence(b2);
		int b3ValenceBefore = MeshOptimization_Apt_Point_Valence(b3);
		int b4ValenceBefore = MeshOptimization_Apt_Point_Valence(b4);
		int b1ValenceAfter = b1ValenceBefore - 1;
		int b2ValenceAfter = b2ValenceBefore - 1;
		int b3ValenceAfter = b3ValenceBefore + 1;
		int b4ValenceAfter = b4ValenceBefore + 1;

		int BeforeValence = 0;
		BeforeValence = abs(b1ValenceBefore - b1ValenceTarget) +
			abs(b2ValenceBefore - b2ValenceTarget) +
			abs(b3ValenceBefore - b3ValenceTarget) +
			abs(b4ValenceBefore - b4ValenceTarget);

		int AfterValence = 0;
		AfterValence = abs(b1ValenceAfter - b1ValenceTarget) +
			abs(b2ValenceAfter - b2ValenceTarget) +
			abs(b3ValenceAfter - b3ValenceTarget) +
			abs(b4ValenceAfter - b4ValenceTarget);
		//object Valence	

		if (AfterValence < BeforeValence) {
			return true;
		}
		else {
			return false;
		}
	}

	void MeshOptimization_Apt_Flip_Processing(int b1, int b2, int b3, int b4) {

		vector<int> b1N = pointNeighbor[b1];
		vector<int> b2N = pointNeighbor[b2];
		vector<int> b3N = pointNeighbor[b3];
		vector<int> b4N = pointNeighbor[b4];

		//update b1 structure
		vector<int> b1N_New;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 == b2 && (b1N3 == b3 || b1N3 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N3 == b3) {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);
				}
				else {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);
				}
			}
			if (b1N3 == b2 && (b1N2 == b3 || b1N2 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N2 == b3) {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);
				}
				else {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);
				}
			}
		}
		if (b1N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b1N_final;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 != -1 && b1N3 != -1) {
				b1N_final.push_back(b1N2);
				b1N_final.push_back(b1N3);
			}
		}
		b1N_final.push_back(b1N_New[0]);
		b1N_final.push_back(b1N_New[1]);

		//update b2 structure
		vector<int> b2N_New;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 == b1 && (b2N3 == b3 || b2N3 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N3 == b3) {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
				else {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
			}
			if (b2N3 == b1 && (b2N2 == b3 || b2N2 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N2 == b3) {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
				else {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
			}
		}
		if (b2N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b2N_final;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 != -1 && b2N3 != -1) {
				b2N_final.push_back(b2N2);
				b2N_final.push_back(b2N3);
			}
		}
		b2N_final.push_back(b2N_New[0]);
		b2N_final.push_back(b2N_New[1]);

		//update b3 structure	
		vector<int> b3N_New;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if ((b3N2 == b1 || b3N2 == b2) && (b3N3 == b1 || b3N3 == b2)) {
				b3N[2 * i] = -1;
				b3N[2 * i + 1] = -1;
				b3N_New.push_back(b3N2);
				b3N_New.push_back(b4);
				b3N_New.push_back(b4);
				b3N_New.push_back(b3N3);
			}
		}
		if (b3N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b3N_final;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if (b3N2 != -1 && b3N3 != -1) {
				b3N_final.push_back(b3N2);
				b3N_final.push_back(b3N3);
			}
		}
		b3N_final.push_back(b3N_New[0]);
		b3N_final.push_back(b3N_New[1]);
		b3N_final.push_back(b3N_New[2]);
		b3N_final.push_back(b3N_New[3]);

		//update b4 structure	
		vector<int> b4N_New;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if ((b4N2 == b1 || b4N2 == b2) && (b4N3 == b1 || b4N3 == b2)) {
				b4N[2 * i] = -1;
				b4N[2 * i + 1] = -1;
				b4N_New.push_back(b4N2);
				b4N_New.push_back(b3);
				b4N_New.push_back(b3);
				b4N_New.push_back(b4N3);
			}
		}
		if (b4N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b4N_final;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if (b4N2 != -1 && b4N3 != -1) {
				b4N_final.push_back(b4N2);
				b4N_final.push_back(b4N3);
			}
		}
		b4N_final.push_back(b4N_New[0]);
		b4N_final.push_back(b4N_New[1]);
		b4N_final.push_back(b4N_New[2]);
		b4N_final.push_back(b4N_New[3]);

		pointNeighbor[b1].clear();
		pointNeighbor[b2].clear();
		pointNeighbor[b3].clear();
		pointNeighbor[b4].clear();

		pointNeighbor[b1] = b1N_final;
		pointNeighbor[b2] = b2N_final;
		pointNeighbor[b3] = b3N_final;
		pointNeighbor[b4] = b4N_final;

	}

	//Relocate vertices on the surface by tangential smoothing
	void MeshOptimization_Apt_TangentSmoothing() {

		double landa = 0.6;
		vector<vector<double>> pointNew = points;

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			vector<double> pi = points[i];
			int p1 = i;
			if (MeshOptimization_Apt_Edge_Point(p1)) {
				pointNew[i] = pi;
				continue;
			}
			vector<int> pointsN = pointNeighbor[i];
			vector<int> pointsN_withoutzRepeat;//neighbor points without repeat point
			vector<double> pointN_Area;//correspoinding area of neighbor point
			vector<double> pointB_List;//Adaptive L ditance for neighbor
			double AreaSum = 0;
			for (int j = 0; j < pointsN.size() / 2; j++) {
				int p2 = pointsN[2 * j];
				int p3 = pointsN[2 * j + 1];
				double areai = MeshOptimization_Apt_TrangularArea(p1, p2, p3);
				//double areai2 = areai / (double)2;
				int indexp2 = MeshOptimization_Apt_Exist_Index(p2, pointsN_withoutzRepeat);
				int indexp3 = MeshOptimization_Apt_Exist_Index(p3, pointsN_withoutzRepeat);
				if (indexp2 == -1) {
					pointsN_withoutzRepeat.push_back(p2);
					pointN_Area.push_back(areai);
				}
				else {
					pointN_Area[indexp2] = (pointN_Area[indexp2] + areai) / 2;
				}
				if (indexp3 == -1) {
					pointsN_withoutzRepeat.push_back(p3);
					pointN_Area.push_back(areai);
				}
				else {
					pointN_Area[indexp3] = (pointN_Area[indexp3] + areai) / 2;
				}
			}


			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {

				double L_bi = MeshOptimization_Apt_KdTreeSearch(pointsN_withoutzRepeat[j]);
				pointB_List.push_back(L_bi);
			}

			for (int j = 0; j < pointN_Area.size(); j++) {

				AreaSum = AreaSum + pointN_Area[j] * pointB_List[j];

			}

			vector<double> gi(3, 0);
			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {
				gi[0] = gi[0] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][0];
				gi[1] = gi[1] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][1];
				gi[2] = gi[2] + ((pointN_Area[j] * pointB_List[j]) / AreaSum) * points[pointsN_withoutzRepeat[j]][2];
			}

			//achieve normal of p1
			vector<double> ni = MeshOptimization_Apt_Normal_Point(p1);
			//achieve move vertor v_pg for gi and pi
			vector<double> v_pg(3);
			v_pg[0] = gi[0] - pi[0];
			v_pg[1] = gi[1] - pi[1];
			v_pg[2] = gi[2] - pi[2];
			double mapN_dis = v_pg[0] * ni[0] + v_pg[1] * ni[1] + v_pg[2] * ni[2];//v_pg map into n distance.
			vector<double> pi_new(3, 0);
			pi_new[0] = pi[0] + landa * (v_pg[0] - mapN_dis * ni[0]);
			pi_new[1] = pi[1] + landa * (v_pg[1] - mapN_dis * ni[1]);
			pi_new[2] = pi[2] + landa * (v_pg[2] - mapN_dis * ni[2]);
			if (!(pi_new[0]<9999 && pi_new[0]>-9999)) {
				continue;
			}
			pointNew[i] = pi_new;
		}

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			points[i][0] = pointNew[i][0];
			points[i][1] = pointNew[i][1];
			points[i][2] = pointNew[i][2];
		}
	}

	//achieve the trangualr area
	double MeshOptimization_Apt_TrangularArea(int b1, int b2, int b3) {
		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> v3(3);
		v3[0] = p1[0] - p3[0];
		v3[1] = p1[1] - p3[1];
		v3[2] = p1[2] - p3[2];

		double v1d = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double v2d = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double v3d = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);

		double p = (v1d + v2d + v3d) / 2;
		double S = sqrt(p * (p - v1d) * (p - v2d) * (p - v3d));

		return S;
	}

	//normal return
	vector<double> MeshOptimization_Apt_Normal_Point(int b1) {

		vector<int> p1N = pointNeighbor[b1];
		vector<double> p1N_Weight;//weights of the normal list
		vector<vector<double>> p1N_Normal;//normal list
		double Area_sum = 0;
		for (int i = 0; i < p1N.size() / 2; i++) {
			int b2 = p1N[2 * i];
			int b3 = p1N[2 * i + 1];
			vector<double> n_i = MeshOptimization_Apt_Normal_Face(b1, b2, b3);
			double area_i = MeshOptimization_Apt_TrangularArea(b1, b2, b3);
			Area_sum = Area_sum + area_i;
			p1N_Weight.push_back(area_i);
			p1N_Normal.push_back(n_i);
		}
		vector<double> n(3, 0);
		for (int i = 0; i < p1N_Weight.size(); i++) {
			n[0] = n[0] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][0];
			n[1] = n[1] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][1];
			n[2] = n[2] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][2];
		}
		return n;
	}

	vector<double> MeshOptimization_Apt_Normal_Face(int b1, int b2, int b3) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3, 0);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3, 0);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> n(3, 0);
		n[0] = v2[2] * v1[1] - v2[1] * v1[2];
		n[1] = -v2[2] * v1[0] + v2[0] * v1[2];
		n[2] = v2[1] * v1[0] - v2[0] * v1[1];

		vector<double> n_Unit = MeshOptimization_Apt_Unit_Normal(n);
		return n_Unit;

	}

	vector<double> MeshOptimization_Apt_Unit_Normal(vector<double> n) {

		double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		n[0] = n[0] / length;
		n[1] = n[1] / length;
		n[2] = n[2] / length;
		return n;

	}

	//Check structure and remove the repeat item from the structure
	vector<int> MeshOptimization_Apt_Structure_RemoveRepeat(vector<int> v) {

		vector<int> result;
		for (int i = 0; i < v.size() / 2; i++) {
			int b1 = v[2 * i];
			int b2 = v[2 * i + 1];
			if (b1 == -1 || b2 == -1) {
				continue;
			}
			for (int j = i + 1; j < v.size() / 2; j++) {
				int b11 = v[2 * j];
				int b12 = v[2 * j + 1];
				if (b11 == -1 || b12 == -1) {
					continue;
				}
				if ((b1 == b11 || b1 == b12) && (b2 == b11 || b2 == b12)) {
					v[2 * j] = -1;
					v[2 * j + 1] = -1;
				}
			}
		}

		for (int i = 0; i < v.size() / 2; i++) {
			int b1 = v[2 * i];
			int b2 = v[2 * i + 1];
			if (b1 == -1 || b2 == -1) {
				continue;
			}
			else {
				result.push_back(b1);
				result.push_back(b2);
			}
		}
		return result;
	}

	//Compute adaptive distance	
	void MeshOptimization_Apt_L_init() {

		cu_ave.clear();
		MeshGeometric mg;
		mg.MeshGeometric_init(points, pointNeighbor);
		cu_ave = mg.MeshGeometric_Get_N();	

		//average
		int step = 3;
		while (step--) {
			vector<double> pL_ave_store = cu_ave;
			for (int i = 0; i < points.size(); i++) {
				vector<int> pointsN = pointNeighbor[i];
				vector<int> pointsNRerepeat;
				vector<double> weight_list;
				double weight_Sum = 0;
				for (int j = 0; j < pointsN.size(); j++) {
					int indexij = pointsN[j];
					if (MeshOptimization_Apt_Exist(indexij, pointsNRerepeat)) {
						continue;
					}
					else {
						pointsNRerepeat.push_back(indexij);
						double weightij = MeshOptimization_Apt_EdgeLength(i, indexij);
						weight_Sum = weight_Sum + weightij;
						weight_list.push_back(weightij);
					}
				}
				double L_ave_i = 0;
				for (int j = 0; j < pointsNRerepeat.size(); j++) {
					L_ave_i = L_ave_i + (weight_list[j] / weight_Sum) * pL_ave_store[pointsNRerepeat[j]];
				}
				cu_ave[i] = L_ave_i;
			}
		}

		double disMax = -9999;
		double disMin = 9999;
		for (int i = 0; i < cu_ave.size(); i++) {
			if (cu_ave[i] < disMin) {
				disMin = cu_ave[i];
			}
			if (cu_ave[i] > disMax) {
				disMax = cu_ave[i];
			}
		}

		

		vector<double> cu_ave_t = cu_ave;
		sort(cu_ave_t.begin(), cu_ave_t.end());
		double li1 = cu_ave_t[cu_ave_t.size() / 5];
		double li2 = cu_ave_t[cu_ave_t.size() * 2 / 5];
		double li3 = cu_ave_t[cu_ave_t.size() * 3 / 5];
		double li4 = cu_ave_t[cu_ave_t.size() * 4 / 5];
		
		double li1_3 = cu_ave_t[cu_ave_t.size() / 3];
		double li2_3 = cu_ave_t[cu_ave_t.size() * 2 / 3];

		/*
		for (int i = 0; i < cu_ave.size(); i++) {

			double rate = (cu_ave[i] - disMin) / (disMax - disMin);
			cu_ave[i] = (rate * 0.4 + 0.8) * L_regular;			

		}*/

		
		for (int i = 0; i < cu_ave.size(); i++) {	
			if (cu_ave[i] < li1) {
				cu_ave[i] = 2 * L_regular;
			}
			else if (cu_ave[i] < li2 && cu_ave[i] >= li1) {
				cu_ave[i] = 1.5 * L_regular;
			}
			else if (cu_ave[i] < li3 && cu_ave[i] >= li2) {
				cu_ave[i] = 1 * L_regular;
			}
			else if (cu_ave[i] < li4 && cu_ave[i] >= li3) {
				cu_ave[i] = 0.8 * L_regular;
			}
			else {
				cu_ave[i] = 0.6 * L_regular;
			}
		}
		

		/*
		for (int i = 0; i < cu_ave.size(); i++) {
			if (cu_ave[i] < li1_3) {
				cu_ave[i] = 1.2 * L_regular;
			}
			else if (cu_ave[i] < li2_3 && cu_ave[i] >= li1_3) {
				cu_ave[i] = 1 * L_regular;
			}			
			else {
				cu_ave[i] = 0.8 * L_regular;
			}
		}
		*/
		

	}
	
	vector<double> MeshOptimization_Apt_Plane(vector<double> p, vector<double> pn,
		vector<double> sourceP) {//map sourceP to the plane

		double A = pn[0];
		double B = pn[1];
		double C = pn[2];
		double D = -A * p[0] - B * p[1] - C * p[2];

		double Q = A * sourceP[0] + B * sourceP[1] + C * sourceP[2];

		double landa = (Q + D) / (A * A + B * B + C * C);

		double xmap = sourceP[0] - A * landa;
		double ymap = sourceP[1] - B * landa;
		double zmap = sourceP[2] - C * landa;

		vector<double> targetP(3);
		targetP[0] = xmap;
		targetP[1] = ymap;
		targetP[2] = zmap;

		return targetP;
	}

	double MeshOptimization_Apt_InnerAngle(int p1, int p2, int p3) {

		vector<double> p1p = points[p1];
		vector<double> p2p = points[p2];
		vector<double> p3p = points[p3];

		vector<double> pn1(3);
		pn1[0] = p2p[0] - p1p[0];
		pn1[1] = p2p[1] - p1p[1];
		pn1[2] = p2p[2] - p1p[2];

		vector<double> pn2(3);
		pn2[0] = p3p[0] - p1p[0];
		pn2[1] = p3p[1] - p1p[1];
		pn2[2] = p3p[2] - p1p[2];

		double cosa1 = pn1[0] * pn2[0] + pn1[1] * pn2[1] + pn1[2] * pn2[2];
		double anglea1 = acos(cosa1);
		return anglea1;


	}

};




