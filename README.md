# multi_stereo_camera
Calibration and localization algorithms for a system of multiple rigidly coupled stereo cameras. All algorithms have been tested on DJI guidance sensor https://www.dji.com/it/guidance. See doc for experimental results.
## Calibration
Calibration algorithm for multi-stereo sensor:
* Extension of http://www-oldurls.inf.ethz.ch/personal/pomarc/pubs/LiIROS13a.pdf
* Based on affine invariant version of SURF
## Localization
Localization based on sparse features reconstruction. Expression of a reconstruction in principal camera frame, localization wrt the reconstruction with three approaches:
* 3D-3D
* Decoupled 3D-2D: PnP over each camera and averaged mean of results
* Coupled 3D-2D: OpenGV (https://laurentkneip.github.io/opengv/) implementation of gPnP
## Visual odometry example
Implementation of a simple VO pipeline (pose only optimization local bundle adjustment).

![Alt Text](https://github.com/emilianogagliardi/multi_stereo_camera/blob/master/doc/demo.gif)
