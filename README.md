## Intrinsic and Isotropic Resampling for 3D Point Clouds

![image](https://user-images.githubusercontent.com/65271555/176191299-220e20ff-a146-48ea-b565-926963d4636c.png)

### Introduction

With rapid development of 3D scanning technology, 3D point cloud based research and applications are becoming more popular. However, major difficulties are still exist which affect the performance of point cloud utilization. Such difficulties include lack of local adjacency information, non-uniform point density, and control of point numbers. In this paper, we propose a two-step intrinsic and isotropic (I\&I) resampling framework to address the challenge of these three major difficulties. The efficient intrinsic control provides geodesic measurement for a point cloud to improve local region detection and avoids redundant geodesic calculation. Then the geometrically-optimized resampling uses a geometric update process to optimize a point cloud into an isotropic or adaptively-isotropic one. The point cloud density can be adjusted to global uniform (isotropic) or local uniform with geometric feature keeping (being adaptively isotropic). The point cloud number can be controlled based on application requirement or user-specification. Experiments show that our point cloud resampling framework achieves outstanding performance in different applications: point cloud simplification, mesh reconstruction and shape registration. 

### Citation
If you find our work useful in your research, please consider citing:

     @article{lv2022intrinsic,
         title={Intrinsic and Isotropic Resampling for 3D Point Clouds},
         author={Lv, Chenlei and Lin, Weisi and Zhao, Baoquan},
         journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
         year={2022},
         publisher={IEEE}
     }

### Update

2022.06.23.

The related paper “Intrinsic and Isotropic Resampling for 3D Point Clouds” has been accepted in TPAMI. 

2021.08.19, Version 0.2.

We upload the program into EXE folder which has been compiled. You can use .exe to test our resampling method directly without any complex configurations. 

2021.01.21, Version 0.1.

The I&I resampling framework for point cloud resampling and mesh reconstruction.

For this version, the geodesic coodinate system is not added into the framework. We use a AIVS method (another work to provide an appximate intrinsic control in simplification) in this version. In futher work, we will upload new codes with full functions. 

The related libiary include:

PCL 1.8.1:  https://pointclouds.org/downloads/

CGAL:       https://www.cgal.org/

freeglut64: http://freeglut.sourceforge.net/index.php#download  
