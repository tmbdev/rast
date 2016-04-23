/**
 * @file   EdgeDetector.cpp 
 * @brief  EdgeDetector in-sight code 
 * @author Son Le (Shawn)
 * @date   15/10/2015
 */

#ifndef __EDGE_DETECTOR_HPP__
#define __EDGE_DETECTOR_HPP__

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <list>
#include <boost/graph/graph_concepts.hpp>

#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/calib3d/calib3d.hpp"



class EdgeDetector
{
public:
  
  struct EdgeData
  {
    cv::Mat                  DxMat;
    cv::Mat                  DyMat;
    cv::Mat                  magMat;
    std::vector<double>      DxArr;
    std::vector<double>      DyArr;
    std::vector<double>      magArr;
    
    /// edgePtsOrg refer to image top-left as origin   
    std::vector<cv::Point2i> edgePtsOrg;
    
    /// edgePts are shifted by COG: coords[m] = temp-COG;
    std::vector<cv::Point2i> edgePts;
  };
  
  EdgeDetector();
  
  ~EdgeDetector();
  
  int EdgeDetect(const cv::Mat image, EdgeData& edDat, cv::Mat& edge_ou);
  
  
private:
  
  cv::Mat Sx, Sy;
  cv::Mat mag, ang;
  cv::Mat img, imgGray;
  
};
#endif



