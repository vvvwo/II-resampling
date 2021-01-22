/*********************************************************************************

						   AIVS-based Remeshing Pro Vision

						      Updating in 2020/09/29

						        By Dr. Chenlei Lv

			The functions includes:
			1. Point cloud Pre-processing;
			2. Isotropic sampling from a point cloud
			3. Anisotropic simplification from a point cloud
			   (Curvature sensitive sampling ot shape feature keeping)
			4. Mesh cropping
			   4.1 remove the trangulars that are corssing the empty box region
			   4.2 adaptive isotropic remeshing
			5. Save the mesh

			Function:
			AIVS_RNPro_PointCloud_Preprocessing(string fileName);
			//point cloud pre-processing

**********************************************************************************/
#pragma once
#include <iostream>
#include "LoadPointCloud.hpp"
#include "normalCompute.hpp"
#include "ballRegionCompute.hpp"
#include "Method_AIVS_SimPro.hpp"
#include "Method_CGAL_Advancing.hpp"
#include "Method_CGAL_Advancing_Edge.hpp"
#include "Mesh_Optimization.hpp"
#include "Mesh_Optimization_Edge.hpp"
#include "Mesh_Optimization_Adaptive.hpp"
#include "Mesh_Optimization_Apt_C.hpp"
#include "Mesh_Optimization_Apt.hpp"
#include "Method_Octree.hpp"
#include <pcl/io/pcd_io.h>  //File input$output
#include <pcl/octree/octree_search.h>  //octree define
#include <pcl/point_types.h>  //point type
using namespace std;
typedef std::vector< pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ> > AlignedPointTVector;
#pragma once

class AIVS_RNPro {

private:

	AIVS_Simplification_Pro asp;
	BallRegion br;
	string fileNameInput;
	double h;
	double weightPointUpdata = 0.6;//update weight for seedpoints

public:

	void AIVS_RNPro_init(string fileName) {

		cout << "pointProcessPipeline_init run!" << endl;
		AIVS_RNPro_PointCloud_Preprocessing(fileName);
		int index = fileName.find_last_of(".");
		fileNameInput = fileName.substr(0, index);		
	
	}

	void AIVS_RNPro_Remesh(string fileNameModel, int pointNum, int classificationIndex, double multi, int optIter) {
		//classificationIndex 1: isotropic; 2: Curvature; 3: Edge.

		cout << "AIVS_Pro_Reconstruct start." << endl;
		vector<vector<double>> seedPointsInit;
		vector<vector<double>> pointCloudPre;

		//1. simplification with upsampling
		if (br.pointCloudData.size() < 3 * pointNum) {
			cout << "0. AIVS_Reconstruct up-sampling:" << endl;
			vector<vector<vector<double>>> result = AIVS_Reconstruct_UpPointsTran(br.pointCloudData, 2);
			cout << "Up-sampling number:" << result[0].size() << endl;
			vector<vector<vector<double>>> result2 = AIVS_RNPro_Octree_DownSampling(result[0], result[1], br.radius / 8);
			pointCloudPre = result2[0];
			cout << "Octree down-sampling number:" << result2.size() << endl;
			vector<int> borderIndex = AIVS_Reconstruct_BorderReturn(result2[0]);
			BallRegion brUp;
			brUp.BallRegion_init(result2[0], result2[1], borderIndex);
			cout << "1. AIVS_Reconstruct Simplification:" << endl;			
			asp.AIVS_Pro_init(brUp, fileNameInput);			
		}
		else if (classificationIndex == 3) {
			
			cout << "1. AIVS_Reconstruct Simplification:" << endl;
			asp.AIVS_Pro_init_Edge(br, fileNameInput);
		
		}
		else {
			cout << "1. AIVS_Reconstruct Simplification:" << endl;
			pointCloudPre = br.pointCloudData;
			asp.AIVS_Pro_init(br, fileNameInput);			
		}

		//classification for different kinds of simplification
		if (classificationIndex == 1) {
			seedPointsInit = asp.AIVS_simplification(pointNum);		
		}
		else if (classificationIndex == 2) {
			vector<double> prate;
			double a1 = 2; double a2 = 3; double a3 = 4; double a4 = 5; double a5 = 6;
			prate.push_back(a1); prate.push_back(a2); prate.push_back(a3);
			prate.push_back(a4); prate.push_back(a5);
			seedPointsInit = asp.AIVS_simplification_Curvature(pointNum, prate);		
		}
		else if (classificationIndex == 3) {
			vector<double> prate;
			prate.push_back(3);
			prate.push_back(7);
			seedPointsInit = asp.AIVS_simplification_Edge_Curvature(pointNum, 3, 7);
			//seedPointsInit = asp.AIVS_simplification_Edge(pointNum, prate);
		}
		else {
			seedPointsInit = asp.AIVS_simplification(pointNum);		
		}

		//2. Init Meshing
		cout << "2. AIVS_Reconstruct Meshing:" << endl;
		vector<vector<int>> faceInforN;
		if (classificationIndex == -1) {
			CGAL_Advancing_R_Edge car;
			car.CGAL_Advancing_R_Edge_init(seedPointsInit);
			car.CGAL_Advancing_R_Edge_FaceInfor();
			seedPointsInit.clear();
			seedPointsInit = car.pointsNewR;
			faceInforN = car.faceInforsR;		
		}
		else {
			CGAL_Advancing_Reconstruct cr;
			cr.CGAL_Advancing_Remesh_init(seedPointsInit);
			faceInforN = cr.CGAL_Advancing_Remesh_FaceInfor();
		}	
		
		
		//4. Check the wrong trangular		
		//cout << "4. AIVS_Reconstruct Face Detect:" << endl;
		//vector<vector<int>> faceInforNT = AIVS_Reconstruct_FaceRemove(faceInforN, seedPointsInit, multi);
		//faceInforN.clear();
		//faceInforN = faceInforNT;
		
		//5. Mesh optimization	
		cout << "5. AIVS_Reconstruct Mesh optimization:" << endl;	
		vector<vector<double>> pointOpt;
		vector<vector<int>> faceOpt;
		if (classificationIndex == 2) {			
			//MeshOptimization_Apt_C moac;
			//moac.MeshOptimization_Apt_C_init(fileNameModel, seedPointsInit, faceInforN);			
			//moac.MeshOptimization_Apt_C_Start(optIter, true);
			//pointOpt = moac.MeshOptimization_Apt_C_GetPoints();
			//faceOpt = moac.MeshOptimization_Apt_C_GetFaces();

			MeshOptimization_Apt mo;
		    mo.MeshOptimization_Apt_init(fileNameModel, seedPointsInit, faceInforN);
			mo.MeshOptimization_Apt_Start(optIter, true);
			pointOpt = mo.MeshOptimization_Apt_GetPoints();
			faceOpt = mo.MeshOptimization_Apt_GetFaces();


		}		
		else if (classificationIndex == 3) {
			MeshOptimization_Edge moe;
			moe.MeshOptimization_init(pointNum, fileNameModel, seedPointsInit, faceInforN);
			moe.MeshOptimization_Start(optIter, true);
			pointOpt = moe.MeshOptimization_GetPoints();
			faceOpt = moe.MeshOptimization_GetFaces();
		}
		else {
			MeshOptimization mo;
			mo.MeshOptimization_init(fileNameModel, seedPointsInit, faceInforN);
			mo.MeshOptimization_Start(optIter, true);
			pointOpt = mo.MeshOptimization_GetPoints();
			faceOpt = mo.MeshOptimization_GetFaces();
		}

		/*
		//5.1 update point position
		for (int i = 0; i < pointOpt.size(); i++) {
			vector<double> pnew = AIVS_pointUpdate(pointOpt[i]);
			pointOpt[i][0] = pnew[0];
			pointOpt[i][1] = pnew[1];
			pointOpt[i][2] = pnew[2];		
		}

		//5.2 search neighbor position
		vector<bool> pointJudSim(br.pointCloudData.size(), false);
		int K = 3;
		for (int i = 0; i < pointOpt.size(); i++) {
			
			//int K = br.pointNumEsti + 1;
			std::vector<int> pointIdxNKNSearch(K);
			std::vector<float> pointNKNSquaredDistance(K);
			std::vector<int> pointNeior;
			//double r_i = pointNKNSquaredDistance[pointNKNSquaredDistance.size() - 1];
			pcl::PointXYZ searchPoint;
			searchPoint.x = pointOpt[i][0];
			searchPoint.y = pointOpt[i][1];
			searchPoint.z = pointOpt[i][2];
			br.kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
			for (int j = 0; j < pointIdxNKNSearch.size(); j++) {
				int indexij = pointIdxNKNSearch[j];
				if (!pointJudSim[indexij]) {
					pointJudSim[indexij] = true;
					break;				
				}			
			}
		}

		pointOpt.clear();
		for (int i = 0; i < pointJudSim.size(); i++) {
			if (pointJudSim[i]) {
				pointOpt.push_back(br.pointCloudData[i]);
			}		
		}*/

		//6. Fix Model
		CGAL_Advancing_Reconstruct Fx;
		Fx.CGAL_Advancing_Remesh_init(pointOpt);
		vector<vector<int>> faceInforNFix = Fx.CGAL_Advancing_Remesh_FaceInfor();

		//7. Store FIle	
		cout << "6. AIVS_Reconstruct Store File:" << endl;		
		string objName_init = "Remesh\\Reconstruction\\" + fileNameModel + to_string(pointNum)+"_"+ "init" + ".obj";	
		string objName_opt = "Remesh\\Reconstruction\\" + fileNameModel + to_string(pointNum) + "_" + "opt" + ".obj";
		string objName_optFix = "Remesh\\Reconstruction\\" + fileNameModel + to_string(pointNum) + "_" + "optFix" + ".obj";
		
		string objName_sim = "Remesh\\Reconstruction\\" + fileNameModel + ".obj";
		//AIVS_RNPro_SaveOBJ(objName_init, seedPointsInit, faceInforN);
		//AIVS_RNPro_SaveOBJ(objName_opt, pointOpt, faceOpt);
		//AIVS_RNPro_SaveOBJ(objName_optFix, pointOpt, faceInforNFix);
		AIVS_RNPro_SaveOBJ(objName_sim, pointOpt, faceInforNFix);
		cout << "AIVS_Reconstruct finish." << endl;	
	}

	void AIVS_RNPro_Remesh_Octree(string fileNameModel, int pointNum, int classificationIndex, double multi) {
		//classificationIndex 1: isotropic; 2: Curvature; 3: Edge.

		cout << "AIVS_Pro_Reconstruct start." << endl;
		vector<vector<double>> seedPointsInit;
		vector<vector<double>> pointCloudPre;

		//1. simplification with upsampling
		
		cout << "0. AIVS_Reconstruct up-sampling:" << endl;
		vector<vector<vector<double>>> result = AIVS_Reconstruct_UpPointsTran(br.pointCloudData, 2);			
		vector<vector<vector<double>>> result2 = AIVS_RNPro_Octree_DownSampling(result[0], result[1], br.radius / 2);
			
		
		
		CGAL_Advancing_Reconstruct cr;
		cr.CGAL_Advancing_Remesh_init(result[0]);
		vector<vector<int>> faceInforN_Upsampling = cr.CGAL_Advancing_Remesh_FaceInfor();

		CGAL_Advancing_Reconstruct cr2;
		cr2.CGAL_Advancing_Remesh_init(result2[0]);
		vector<vector<int>> faceInforN_Octree = cr2.CGAL_Advancing_Remesh_FaceInfor();
		

		//6. Store FIle	
		cout << "6. AIVS_Reconstruct Store File:" << endl;
		string objName1 = "Octree\\" + fileNameModel + "_1.obj";
		string objName2 = "Octree\\" + fileNameModel + "_2.obj";
		

		//AIVS_SaveOFF(offName, seedPoints, faceInforN);
		//if (classificationIndex == 3) {
			//AIVS_RNPro_SaveOBJ(objNameE, seedPointsInit, faceInforEdge);		
		//}

		AIVS_RNPro_SaveOBJ(objName1, result[0], faceInforN_Upsampling);
		AIVS_RNPro_SaveOBJ(objName2, result2[0], faceInforN_Octree);
		//AIVS_RNPro_SavePoints(objNameP, pointCloudPre);
		//AIVS_RNPro_SavePoints(objNameP, pointEdges);
		cout << "AIVS_Reconstruct finish." << endl;
	}

private:

	vector<vector<vector<double>>> AIVS_Reconstruct_UpPointsTran(vector<vector<double>> p, int Up) {

		cout << "increase the points number from the point cloud" << endl;

		//1. achieve mesh
		CGAL_Advancing_Reconstruct crt;
		crt.CGAL_Advancing_Remesh_init(p);
		vector<vector<int>> faceInfor = crt.CGAL_Advancing_Remesh_FaceInfor();

		//2. achieve regular area
		vector<int> faceInfor1 = faceInfor[1];
		vector<vector<double>> b1;
		b1.push_back(p[faceInfor1[0]]);
		b1.push_back(p[faceInfor1[1]]);
		b1.push_back(p[faceInfor1[2]]);

		vector<int> faceInfor2 = faceInfor[2];
		vector<vector<double>> b2;
		b2.push_back(p[faceInfor2[0]]);
		b2.push_back(p[faceInfor2[1]]);
		b2.push_back(p[faceInfor2[2]]);

		vector<int> faceInfor3 = faceInfor[3];
		vector<vector<double>> b3;
		b3.push_back(p[faceInfor3[0]]);
		b3.push_back(p[faceInfor3[1]]);
		b3.push_back(p[faceInfor3[2]]);

		double avergeN = (AIVS_Reconstruct_AchieveArea(b1) +
			AIVS_Reconstruct_AchieveArea(b2) + AIVS_Reconstruct_AchieveArea(b3)) / 3;
		double points = (1 + Up) * Up / 2;

		//3. achieve regular area
		for (int i = 0; i < faceInfor.size(); i++) {
			vector<int> faceInfor1 = faceInfor[i];
			vector<vector<double>> bi;
			bi.push_back(p[faceInfor1[0]]);
			bi.push_back(p[faceInfor1[1]]);
			bi.push_back(p[faceInfor1[2]]);
			double areai = AIVS_Reconstruct_AchieveArea(bi);
			int N = points * areai / avergeN;
			vector<vector<double>> pnewPoint = AIVS_Reconstruct_NewPoint(bi, N);
			p.insert(p.end(), pnewPoint.begin(), pnewPoint.end());
		}

		//4. update points positions 
		vector<vector<double>> pNormal(p.size());

#pragma omp parallel for
		for (int i = 0; i < p.size(); i++) {
			if (i % 1000 == 0) {
				cout << (p.size() - i) / 1000 << " ";
			}
			vector<double> pNewi = AIVS_pointUpdate(p[i]);
			p[i][0] = pNewi[0];
			p[i][1] = pNewi[1];
			p[i][2] = pNewi[2];
			vector<double> pNi(3);
			pNi[0] = pNewi[3];
			pNi[1] = pNewi[4];
			pNi[2] = pNewi[5];
			pNormal[i] = pNi;
		}

		//remove repeat points
		vector<vector<double>> pRemoveRepeat;
		vector<vector<double>> pRemoveRepeatN;

		int kn = 14;
		pcl::KdTreeFLANN<pcl::PointXYZ> kdtreeSeed = AIVS_KdTreeConstruct(p);//kdtree
		for (int i = 0; i < p.size(); i++) {
			if ((p.size() - i) % 1000 == 0) {
				cout << (p.size() - i) / 1000 << " ";
			}
			if (p[i][0] == -9999) {
				continue;
			}
			else {
				vector<double> seedPoints_i = p[i];
				vector<int> pointIdxNKNSearch(kn);
				vector<float> pointNKNSquaredDistance(kn);
				std::vector<int> pointNeibor_i;
				pcl::PointXYZ searchPoint;
				searchPoint.x = seedPoints_i[0];
				searchPoint.y = seedPoints_i[1];
				searchPoint.z = seedPoints_i[2];
				kdtreeSeed.nearestKSearch(searchPoint, kn, pointIdxNKNSearch, pointNKNSquaredDistance);

				if (pointNKNSquaredDistance[1] != 0) {
					pRemoveRepeat.push_back(p[i]);
					pRemoveRepeatN.push_back(pNormal[i]);
				}
				else {
					int index = pointIdxNKNSearch[1];
					if (i > index) {
						//cout << "Error! i = " << i << endl;					
					}
					else {
						pRemoveRepeat.push_back(p[i]);
						pRemoveRepeatN.push_back(pNormal[i]);
						p[index][0] = -9999;
					}
				}
			}
		}
		cout << endl;
		vector<vector<vector<double>>> result;
		vector<vector<double>> keyPoints;
		result.push_back(pRemoveRepeat);
		result.push_back(pRemoveRepeatN);
		result.push_back(keyPoints);

		return result;

	}

	double AIVS_Reconstruct_AchieveArea(vector<vector<double>> vPoint) {

		double x0 = vPoint[0][0];
		double y0 = vPoint[0][1];
		double z0 = vPoint[0][2];

		double x1 = vPoint[1][0];
		double y1 = vPoint[1][1];
		double z1 = vPoint[1][2];

		double x2 = vPoint[2][0];
		double y2 = vPoint[2][1];
		double z2 = vPoint[2][2];

		double b1 = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
		double b2 = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2) + (z0 - z2) * (z0 - z2));
		double b3 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));

		double p = (b1 + b2 + b3) / 2;
		double tArea = sqrt(p * (p - b1) * (p - b2) * (p - b3));
		return tArea;
	}

	vector<vector<double>> AIVS_Reconstruct_NewPoint(vector<vector<double>> vPoint, int N) {
		//N is area parameter
		vector<vector<double>> newPoint;

		//init three points in trangular 
		double x0 = vPoint[0][0];
		double y0 = vPoint[0][1];
		double z0 = vPoint[0][2];

		double x1 = vPoint[1][0];
		double y1 = vPoint[1][1];
		double z1 = vPoint[1][2];

		double x2 = vPoint[2][0];
		double y2 = vPoint[2][1];
		double z2 = vPoint[2][2];

		//compute border point 
		int Q = sqrt(2 * N + 0.5) - 1.414 / 2 + 1;

		vector<vector<double>> b1;
		vector<vector<double>> b2;

		for (int i = 1; i < Q; i++) {

			vector<double> b1i(3);
			b1i[0] = x0 * (1 - (double)i / (double)Q) + x1 * ((double)i / (double)Q);
			b1i[1] = y0 * (1 - (double)i / (double)Q) + y1 * ((double)i / (double)Q);
			b1i[2] = z0 * (1 - (double)i / (double)Q) + z1 * ((double)i / (double)Q);
			b1.push_back(b1i);

			vector<double> b2i(3);
			b2i[0] = x0 * (1 - (double)i / (double)Q) + x2 * ((double)i / (double)Q);
			b2i[1] = y0 * (1 - (double)i / (double)Q) + y2 * ((double)i / (double)Q);
			b2i[2] = z0 * (1 - (double)i / (double)Q) + z2 * ((double)i / (double)Q);
			b2.push_back(b2i);

		}
		vector<vector<double>> pMiddlePoints;
		for (int i = 0; i < b1.size(); i++) {
			int seg = i + 2;
			vector<double> b1i = b1[i];
			vector<double> b2i = b2[i];
			for (int j = 0; j < i + 1; j++) {
				double weightij = (double)(j + 1) / (double)seg;
				vector<double> newij(3);
				newij[0] = b1i[0] * (1 - weightij) + b2i[0] * weightij;
				newij[1] = b1i[1] * (1 - weightij) + b2i[1] * weightij;
				newij[2] = b1i[2] * (1 - weightij) + b2i[2] * weightij;
				newPoint.push_back(newij);
			}
		}
		newPoint.insert(newPoint.end(), b1.begin(), b1.end());
		newPoint.insert(newPoint.end(), b2.begin(), b2.end());
		return newPoint;
	}

	void AIVS_RNPro_PointCloud_Preprocessing(string fileName) {			

		//load program, support obj, off, ply, xyz, txt	
		LoadPointCloud lpc;
		lpc.PointCloud_Load(fileName);
		string fileNameNormal = lpc.FileNormal;
		vector<vector<double>> pdata = lpc.pointSet_uniform;
		//octree down-sampling if needed
		vector<int> borderIndex = lpc.indexBorder;
		if (pdata.size() >= 600000) {
			vector<vector<double>> p_d = AIVS_RNPro_Octree_DownSampling(pdata);
			pdata.clear();
			pdata = p_d;
			borderIndex.clear();
			borderIndex =  AIVS_ReturnBorderIndex(pdata);			
		}
		//normal estimation
		NormalEstimation ne;
		ne.estimateNormal_init(fileNameNormal);
		if (ne.normalLoad()) {
			cout << "pointProcessPipeline_init: normal file exist and the data are loaded." << endl;
		}		
		else {
			cout << "pcl normal start." << endl;				
			ne.estimateNormal_PCL_MP(pdata);
		}		
		vector<vector<double>> ndata = ne.normalVector;		
		br.BallRegion_init(pdata, ndata, borderIndex);	
	}

	vector<vector<double>> AIVS_RNPro_Octree_DownSampling(vector<vector<double>> pData) {
		cout << "Octree down-sampling start:"<<endl;
		PCL_octree po;
		vector<vector<double>> point_down_sampling = po.PCL_Octree_Simplification_WithOutNormal(pData);
		cout << "Octree down-sampling point: " << pData.size() << " to " << point_down_sampling.size() << endl;
		return point_down_sampling;
	}

	vector<vector<vector<double>>> AIVS_RNPro_Octree_DownSampling(vector<vector<double>> pData,
		vector<vector<double>> nData, float resolution) {

		cout << "PCL down-sampling start:" << endl;
		clock_t t1;
		clock_t t2;
		t1 = clock();
		vector<vector<double>> resultP;
		vector<vector<double>> resultN;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		// fills a PointCloud with random data
		for (int i = 0; i < pData.size(); i++)
		{
			pcl::PointXYZ cloud_i;
			cloud_i.x = pData[i][0];
			cloud_i.y = pData[i][1];
			cloud_i.z = pData[i][2];
			cloud->push_back(cloud_i);
		}
		// construct orctree
		//float resolution = 0.03f;  //设置octree体素分辨率
		pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(resolution); //建立octree对象
		octree.setInputCloud(cloud); //传入需要建立kdtree的点云指针
		octree.addPointsFromInputCloud();  //构建Octree		
		AlignedPointTVector voxel_center_list_arg;
		octree.getOccupiedVoxelCenters(voxel_center_list_arg);
		for (int i = 0; i < voxel_center_list_arg.size(); i++) {
			pcl::PointXYZ pdi = voxel_center_list_arg[i];
			int K = 1;
			std::vector<int > pointIdxNKNSearch;
			std::vector<float> pointNKNSquaredDistance;
			octree.nearestKSearch(pdi, K, pointIdxNKNSearch, pointNKNSquaredDistance);
			vector<double> pi;
			vector<double> ni;
			pi.push_back(pData[pointIdxNKNSearch[0]][0]);
			pi.push_back(pData[pointIdxNKNSearch[0]][1]);
			pi.push_back(pData[pointIdxNKNSearch[0]][2]);
			ni.push_back(nData[pointIdxNKNSearch[0]][0]);
			ni.push_back(nData[pointIdxNKNSearch[0]][1]);
			ni.push_back(nData[pointIdxNKNSearch[0]][2]);
			resultP.push_back(pi);
			resultN.push_back(ni);
		}
		t2 = clock();
		std::cout << "Octree down-sampling:" << (t2 - t1) / 1000.0 << "s" << endl;
		cout << "Original Point:" << pData.size() << endl;
		cout << "Down-sampling Point:" << resultP.size() << endl;
		vector<vector<vector<double>>> resultFinal;
		resultFinal.push_back(resultP);
		resultFinal.push_back(resultN);
		return resultFinal;
	}

	vector<vector<double>> AIVS_RNPro_Octree_DownSampling_ChangePosition(vector<vector<double>> pData,
		float resolution) {

		cout << "PCL down-sampling start:" << endl;
		clock_t t1;
		clock_t t2;
		t1 = clock();
		vector<vector<double>> resultP;
		vector<vector<double>> resultN;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		// fills a PointCloud with random data
		for (int i = 0; i < pData.size(); i++)
		{
			pcl::PointXYZ cloud_i;
			cloud_i.x = pData[i][0];
			cloud_i.y = pData[i][1];
			cloud_i.z = pData[i][2];
			cloud->push_back(cloud_i);
		}
		// construct orctree
		//float resolution = 0.03f;  //设置octree体素分辨率
		pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree(resolution); //建立octree对象
		octree.setInputCloud(cloud); //传入需要建立kdtree的点云指针
		octree.addPointsFromInputCloud();  //构建Octree		
		AlignedPointTVector voxel_center_list_arg;
		octree.getOccupiedVoxelCenters(voxel_center_list_arg);
		for (int i = 0; i < voxel_center_list_arg.size(); i++) {
			pcl::PointXYZ pdi = voxel_center_list_arg[i];			
			vector<double> pi;			
			pi.push_back(pdi.x);
			pi.push_back(pdi.y);
			pi.push_back(pdi.z);
			resultP.push_back(pi);			
		}		
		return resultP;
	}

	void AIVS_RNPro_SaveOBJ(string fileName, vector<vector<double>> points, vector<vector<int>> facet) {

		ofstream f1(fileName);

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < facet.size(); i++) {
			f1 << "f " << facet[i][0] + 1 << " " << facet[i][1] + 1 << " " << facet[i][2] + 1 << endl;
		}

		f1.close();


	}

	void AIVS_RNPro_SavePoints(string fileName, vector<vector<double>> points) {

		ofstream f1(fileName);

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}		

		f1.close();
	}

	pcl::KdTreeFLANN<pcl::PointXYZ> AIVS_KdTreeConstruct(vector<vector<double>> seedpoints) {

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

	vector<double> AIVS_pointUpdate(vector<double> point_i) {

		int iter = 10;

#pragma region Achieve Neibor
		//Achieve Neibor
		int K = 13;
		//int K = br.pointNumEsti + 1;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		std::vector<int> pointNeior;
		//double r_i = pointNKNSquaredDistance[pointNKNSquaredDistance.size() - 1];
		pcl::PointXYZ searchPoint;
		searchPoint.x = point_i[0];
		searchPoint.y = point_i[1];
		searchPoint.z = point_i[2];
		br.kdtree.nearestKSearch(searchPoint, K, pointIdxNKNSearch, pointNKNSquaredDistance);
		if (pointNKNSquaredDistance[0] == 0) {
			point_i.push_back(br.pointNormal[pointIdxNKNSearch[0]][0]);
			point_i.push_back(br.pointNormal[pointIdxNKNSearch[0]][1]);
			point_i.push_back(br.pointNormal[pointIdxNKNSearch[0]][2]);
			return point_i;
		}
		else {
			pointNeior.insert(pointNeior.end(), pointIdxNKNSearch.begin(), pointIdxNKNSearch.end());
		}

		vector<vector<double>> pointNormal_i(pointNeior.size());
		for (int i = 0; i < pointIdxNKNSearch.size(); i++) {
			pointNormal_i[i] = br.pointNormal[pointIdxNKNSearch[i]];
		}
		//regular normal
		vector<double> normalUniform = pointNormal_i[0];
		for (int i = 0; i < pointNormal_i.size(); i++) {
			vector<double> normalNeibor_i = pointNormal_i[i];
			vector<double> normalNeibor_i2;
			normalNeibor_i2.push_back(-normalNeibor_i[0]);
			normalNeibor_i2.push_back(-normalNeibor_i[1]);
			normalNeibor_i2.push_back(-normalNeibor_i[2]);
			//manage different directions
			double a1 = normalUniform[0] * normalNeibor_i[0] +
				normalUniform[1] * normalNeibor_i[1] + normalUniform[2] * normalNeibor_i[2];
			double a2 = normalUniform[0] * normalNeibor_i2[0] +
				normalUniform[1] * normalNeibor_i2[1] + normalUniform[2] * normalNeibor_i2[2];
			double cosa1 = acos(a1);
			double cosa2 = acos(a2);
			if (cosa1 > cosa2) {
				pointNormal_i[i][0] = normalNeibor_i2[0];
				pointNormal_i[i][1] = normalNeibor_i2[1];
				pointNormal_i[i][2] = normalNeibor_i2[2];
			}
		}

#pragma endregion

#pragma region MLS error
		//++++++++++++++++++++interater start+++++++++++++++++++++++++++
		double errorExist = 0.0001;
		vector<double> px;
		px.insert(px.end(), point_i.begin(), point_i.end());
		//vector<int> p_neibor = br.pointNeibor[i];			
		//regularNoraml

		vector<double> px_store(3);
		vector<double> nx_store(3);
		vector<double> ax;//new point position
		ax.push_back(0);
		ax.push_back(0);
		ax.push_back(0);
		vector<double> nx;//new point normal
		nx.push_back(0);
		nx.push_back(0);
		nx.push_back(0);
		double errorEndTem;//record new 
		double errorStore = 9999;
		double weight;
		while (iter) {
			//vector<double> ax = simMeasurement_cop_a(px, pointNeiborRegualrNum);
			//vector<double> nx = simMeasurement_cop_n(px, pointNeiborRegualrNum, pointNeiborNormalRegualrNum);
			double fenmu = 0;
			for (int j = 0; j < pointNeior.size(); j++) {
				double dis_i = sqrt((br.pointCloudData[pointNeior[j]][0] - px[0]) * (br.pointCloudData[pointNeior[j]][0] - px[0]) +
					(br.pointCloudData[pointNeior[j]][1] - px[1]) * (br.pointCloudData[pointNeior[j]][1] - px[1]) +
					(br.pointCloudData[pointNeior[j]][2] - px[2]) * (br.pointCloudData[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				double eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + br.pointCloudData[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + br.pointCloudData[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + br.pointCloudData[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointNormal_i[j][0] * eData;
				nx[1] = nx[1] + pointNormal_i[j][1] * eData;
				nx[2] = nx[2] + pointNormal_i[j][2] * eData;
				fenmu = fenmu + eData;
			}
			if (fenmu != 0) {
				ax[0] = ax[0] / fenmu;
				ax[1] = ax[1] / fenmu;
				ax[2] = ax[2] / fenmu;
				nx[0] = nx[0] / fenmu;
				nx[1] = nx[1] / fenmu;
				nx[2] = nx[2] / fenmu;
			}
			fenmu = 0;

			//5.3 Set x' = x - n(x')T(a(x')-x)n(x'), weight = n(x')T(a(x')-x)
			weight = nx[0] * (point_i[0] - ax[0]) +
				nx[1] * (point_i[1] - ax[1]) + nx[2] * (point_i[2] - ax[2]);
			px_store[0] = px[0];
			px_store[1] = px[1];
			px_store[2] = px[2];
			nx_store[0] = nx[0];
			nx_store[1] = nx[1];
			nx_store[2] = nx[2];
			px[0] = point_i[0] - weight * nx[0];
			px[1] = point_i[1] - weight * nx[1];
			px[2] = point_i[2] - weight * nx[2];
			//5.4 ||n(x')T(a(x')-x)n(x')||>errorEndTem
			ax[0] = 0;
			ax[1] = 0;
			ax[2] = 0;
			nx[0] = 0;
			nx[1] = 0;
			nx[2] = 0;

			//vector<double> axnew = simMeasurement_cop_a(px, pointNeiborRegualrNum);
			//vector<double> nxnew = simMeasurement_cop_n(px, pointNeiborRegualrNum, pointNeiborNormalRegualrNum);

			for (int j = 0; j < pointNeior.size(); j++) {
				double dis_i = sqrt((br.pointCloudData[pointNeior[j]][0] - px[0]) * (br.pointCloudData[pointNeior[j]][0] - px[0]) +
					(br.pointCloudData[pointNeior[j]][1] - px[1]) * (br.pointCloudData[pointNeior[j]][1] - px[1]) +
					(br.pointCloudData[pointNeior[j]][2] - px[2]) * (br.pointCloudData[pointNeior[j]][2] - px[2]));
				if (dis_i == 0) {
					continue;
				}
				double eData = -((dis_i / h) * (dis_i / h));
				eData = exp(eData);
				ax[0] = ax[0] + br.pointCloudData[pointNeior[j]][0] * eData;
				ax[1] = ax[1] + br.pointCloudData[pointNeior[j]][1] * eData;
				ax[2] = ax[2] + br.pointCloudData[pointNeior[j]][2] * eData;
				nx[0] = nx[0] + pointNormal_i[j][0] * eData;
				nx[1] = nx[1] + pointNormal_i[j][1] * eData;
				nx[2] = nx[2] + pointNormal_i[j][2] * eData;
				fenmu = fenmu + eData;
			}
			if (fenmu != 0) {
				ax[0] = ax[0] / fenmu;
				ax[1] = ax[1] / fenmu;
				ax[2] = ax[2] / fenmu;
				nx[0] = nx[0] / fenmu;
				nx[1] = nx[1] / fenmu;
				nx[2] = nx[2] / fenmu;
			}

			weight = nx[0] * (point_i[0] - ax[0]) +
				nx[1] * (point_i[1] - ax[1]) + nx[2] * (point_i[2] - ax[2]);

			ax[0] = 0;
			ax[1] = 0;
			ax[2] = 0;
			nx[0] = 0;
			nx[1] = 0;
			nx[2] = 0;

			errorEndTem = abs(weight);
			if (errorEndTem < errorStore) {
				errorStore = errorEndTem;
			}
			else {
				px[0] = px_store[0];
				px[1] = px_store[1];
				px[2] = px_store[2];
				break;
			}
			if (errorEndTem < errorExist) {//|| errorEndTem > errorstore
				break;
			}
			iter--;
		}
#pragma endregion 

		vector<double> finalResult(6);
		finalResult[0] = px[0];
		finalResult[1] = px[1];
		finalResult[2] = px[2];
		finalResult[3] = nx_store[0];
		finalResult[4] = nx_store[1];
		finalResult[5] = nx_store[2];
		return finalResult;
	}

	vector<int> AIVS_Reconstruct_BorderReturn(vector<vector<double>> p) {

		vector<int> indexBorder;

		double maxx, maxy, maxz, minx, miny, minz;

		maxx = minx = p[0][0];
		maxy = miny = p[0][1];
		maxz = minz = p[0][2];

		int indexMaxX = -1;
		int indexMaxY = -1;
		int indexMaxZ = -1;
		int indexMinX = -1;
		int indexMinY = -1;
		int indexMinZ = -1;

		for (int i = 0; i < p.size(); i++) {
			double xi = p[i][0];
			double yi = p[i][1];
			double zi = p[i][2];
			if (xi < minx) {
				minx = xi;
				indexMinX = i;
			}
			if (xi > maxx) {
				maxx = xi;
				indexMaxX = i;
			}
			if (yi < miny) {
				miny = yi;
				indexMinY = i;
			}
			if (yi > maxy) {
				maxy = yi;
				indexMaxY = i;
			}
			if (zi < minz) {
				minz = zi;
				indexMinZ = i;
			}
			if (zi > maxz) {
				maxz = zi;
				indexMaxZ = i;
			}
		}
		indexBorder.push_back(indexMinX);
		indexBorder.push_back(indexMinY);
		indexBorder.push_back(indexMinZ);
		indexBorder.push_back(indexMaxX);
		indexBorder.push_back(indexMaxY);
		indexBorder.push_back(indexMaxZ);
		return indexBorder;
	}

	//Close holes
	vector<vector<int>> AIVS_CloseHoles(vector<vector<int>> neighbor) {

		//inner border
		vector<vector<int>> holesBorders;
		vector<vector<int>> faceBorders;
		for (int i = 0; i < neighbor.size(); i++) {
			int b1 = i;
			vector<int> faceInfor_i = neighbor[i];
			vector<int> faceInfor_i_Repeat;
			vector<int> faceInfor_i_Num;
			vector<int> holes;
			for (int j = 0; j < faceInfor_i.size()/2; j++) {
				int b2 = faceInfor_i[2 * j];
				int b3 = faceInfor_i[2 * j + 1];
				int b2Index = AIVS_ExistVector(b2, faceInfor_i_Repeat);
				int b3Index = AIVS_ExistVector(b3, faceInfor_i_Repeat);
				if (b2Index == -1) {
					faceInfor_i_Repeat.push_back(b2);
					faceInfor_i_Num.push_back(1);
				}
				else {
					faceInfor_i_Num[b2Index]++;
				}
				if (b3Index == -1) {
					faceInfor_i_Repeat.push_back(b3);
					faceInfor_i_Num.push_back(1);
				}
				else {
					faceInfor_i_Num[b3Index]++;				
				}
			}

			for (int j = 0; j < faceInfor_i_Num.size(); j++) {
				int indexNum = faceInfor_i_Num[j];
				if (indexNum < 2) {
					holes.push_back(faceInfor_i_Repeat[j]);
					for (int k = 0; k < faceInfor_i.size()/2; k++) {
						int fk2 = faceInfor_i[2 * k];
						int fk3 = faceInfor_i[2 * k + 1];
						if (fk2 == faceInfor_i_Repeat[j] || fk3 == faceInfor_i_Repeat[j]) {
							vector<int> faceBorder_i(3);
							faceBorder_i[0] = i;
							faceBorder_i[1] = fk2;
							faceBorder_i[2] = fk3;
							faceBorders.push_back(faceBorder_i);
							break;
						}
					}
				}			
			}

			for (int j = 0; j < holes.size(); j++) {				
				if (b1 > holes[j]) {
					continue;				
				}
				else {					
					vector<int> holesBorders_i;
					holesBorders_i.push_back(b1);
					holesBorders_i.push_back(holes[j]);
					holesBorders.push_back(holesBorders_i);
				}			
			}
		}

		vector<bool> holesBordersJudge(holesBorders.size(), false);
		vector<vector<int>> borderPointsList;
		for (int i = 0; i < holesBorders.size(); i++) {
			if (holesBordersJudge[i]) {
				continue;			
			}
			else {
				holesBordersJudge[i] = true;			
			}
			vector<int> borderPoints;
			borderPoints.push_back(holesBorders[i][0]);
			borderPoints.push_back(holesBorders[i][1]);

			int searchIndex = holesBorders[i][1];

			while (1) {
				bool processingJudge = false;
				for (int j = 0; j < holesBorders.size(); j++) {	
					if (holesBordersJudge[j]) {
						continue;					
					}
					int h2j = holesBorders[j][0];
					int h3j = holesBorders[j][1];
					if (h2j == searchIndex && AIVS_ExistVector(h3j, borderPoints) == -1) {
						processingJudge = true;
						holesBordersJudge[j] = true;
						borderPoints.push_back(h3j);
						searchIndex = h3j;
					}
					else if (h3j == searchIndex && AIVS_ExistVector(h2j, borderPoints) == -1) {
						processingJudge = true;
						holesBordersJudge[j] = true;
						borderPoints.push_back(h2j);
						searchIndex = h2j;
					}
					else {
						continue;					
					}
				}
				if (!processingJudge) {
					break;				
				}
			}
			if (borderPoints.size() > 2) {
				borderPointsList.push_back(borderPoints);			
			}			
		}	

		//Rebuild:
		vector<vector<int>> faceNew;
		for (int i = 0; i < borderPointsList.size(); i++) {
			vector<int> borderPoints_i = borderPointsList[i];
			int startPoints = borderPoints_i[0];
			vector<int> faceInfor_i;
			faceInfor_i.push_back(startPoints);
			while (1) {
				bool judge = false;
				for (int j = 0; j < holesBorders.size(); j++) {
					int b2 = holesBorders[j][0];
					int b3 = holesBorders[j][1];
					if (b2 == startPoints && AIVS_ExistVector(b3, faceInfor_i) == -1) {
						judge = true;
						faceInfor_i.push_back(b3);
						startPoints = b3;
					}
					else if(b3 == startPoints && AIVS_ExistVector(b2, faceInfor_i) == -1){
						judge = true;
						faceInfor_i.push_back(b2);
						startPoints = b2;					
					}
					else {
						continue;					
					}				
				}
				if (!judge) {
					break;				
				}
			}
			//startPoints
			for (int j = 1; j < faceInfor_i.size(); j++) {
				int startPoints2 = faceInfor_i[j - 1];
				int startPoints3 = faceInfor_i[j];
				for (int k = 0; k < faceBorders.size(); k++) {
					int k1 = faceBorders[k][0];
					int k2 = faceBorders[k][1];
					int k3 = faceBorders[k][2];
					if (k1 == startPoints2 || k1 == startPoints3) {
						if (k2 == startPoints2 || k2 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k2;
							face_New_i[2] = k1;
							faceNew.push_back(face_New_i);	
							break;
						}
						else if (k3 == startPoints2 || k3 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k1;
							face_New_i[2] = k3;
							faceNew.push_back(face_New_i);
							break;						
						}					
					}
					else if (k2 == startPoints2 || k2== startPoints3) {
						if (k1 == startPoints2 || k1 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k2;
							face_New_i[2] = k1;
							faceNew.push_back(face_New_i);
							break;
						}
						else if (k3 == startPoints2 || k3 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k3;
							face_New_i[2] = k2;
							faceNew.push_back(face_New_i);
							break;

						}					
					}
					else if (k3 == startPoints2 || k3 == startPoints3) {
						if (k1 == startPoints2 || k1 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k1;
							face_New_i[2] = k3;
							faceNew.push_back(face_New_i);
							break;

						}
						else if (k2 == startPoints2 || k2 == startPoints3) {
							vector<int> face_New_i(3);
							face_New_i[0] = startPoints;
							face_New_i[1] = k3;
							face_New_i[2] = k2;
							faceNew.push_back(face_New_i);
							break;
						}
					}				
				}			
			}		
		}

		return faceNew;
	}

	vector<vector<int>> AIVS_Reconstruct_FaceRemove(vector<vector<int>> faceInfor, vector<vector<double>> seedPoints, double multi) {

		double scale = br.unitSize * multi;
		vector<vector<int>> faceInforN;
		for (int i = 0; i < faceInfor.size(); i++) {

			vector<double> b0(3);
			b0[0] = seedPoints[faceInfor[i][0]][0];
			b0[1] = seedPoints[faceInfor[i][0]][1];
			b0[2] = seedPoints[faceInfor[i][0]][2];

			vector<double> b1(3);
			b1[0] = seedPoints[faceInfor[i][1]][0];
			b1[1] = seedPoints[faceInfor[i][1]][1];
			b1[2] = seedPoints[faceInfor[i][1]][2];

			vector<double> b2(3);
			b2[0] = seedPoints[faceInfor[i][2]][0];
			b2[1] = seedPoints[faceInfor[i][2]][1];
			b2[2] = seedPoints[faceInfor[i][2]][2];

			double b01 = sqrt((b0[0] - b1[0]) * (b0[0] - b1[0])
				+ (b0[1] - b1[1]) * (b0[1] - b1[1])
				+ (b0[2] - b1[2]) * (b0[2] - b1[2]));

			double b12 = sqrt((b2[0] - b1[0]) * (b2[0] - b1[0])
				+ (b2[1] - b1[1]) * (b2[1] - b1[1])
				+ (b2[2] - b1[2]) * (b2[2] - b1[2]));

			double b02 = sqrt((b0[0] - b2[0]) * (b0[0] - b2[0])
				+ (b0[1] - b2[1]) * (b0[1] - b2[1])
				+ (b0[2] - b2[2]) * (b0[2] - b2[2]));

			
			if (b01 > scale || b12 > scale || b02 > scale)
			{
				continue;
			}
			else {
				faceInforN.push_back(faceInfor[i]);
			}
		}
		return faceInforN;


	}

	int AIVS_ExistVector(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return i;			
			}		
		}
		return -1;
	
	}

	vector<int> AIVS_ReturnBorderIndex(vector<vector<double>> pointSet) {
		
		vector<int> indexBorder;

		double maxx, maxy, maxz, minx, miny, minz;
		
		maxx = minx = pointSet[0][0];
		maxy = miny = pointSet[0][1];
		maxz = minz = pointSet[0][2];

		int indexMaxX = 0;
		int indexMaxY = 0;
		int indexMaxZ = 0;
		int indexMinX = 0;
		int indexMinY = 0;
		int indexMinZ = 0;

		for (int i = 0; i < pointSet.size(); i++) {
			double xi = pointSet[i][0];
			double yi = pointSet[i][1];
			double zi = pointSet[i][2];
			if (xi < minx) {
				minx = xi;
				indexMinX = i;
			}
			if (xi > maxx) {
				maxx = xi;
				indexMaxX = i;
			}
			if (yi < miny) {
				miny = yi;
				indexMinY = i;
			}
			if (yi > maxy) {
				maxy = yi;
				indexMaxY = i;
			}
			if (zi < minz) {
				minz = zi;
				indexMinZ = i;
			}
			if (zi > maxz) {
				maxz = zi;
				indexMaxZ = i;
			}
		}

		indexBorder.push_back(indexMinX);
		indexBorder.push_back(indexMinY);
		indexBorder.push_back(indexMinZ);
		indexBorder.push_back(indexMaxX);
		indexBorder.push_back(indexMaxY);
		indexBorder.push_back(indexMaxZ);	
		return indexBorder;
	}

};
