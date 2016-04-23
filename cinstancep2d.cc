/* Copyright (c) 1990-1995 by Thomas M. Breuel */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>


#include "misc.h"
#include "narray.h"
#include "vec2.h"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "EdgeDetector.hpp"

using namespace colib;

#include "util.h"
#include "rast.h"

namespace lumo_cinstancep2d {


inline int urand48() { return abs(int(lrand48())); }

#if 0
    template <class T>
    static void shuffle(narray<T> &narray) {
        int n = narray.length();
        for(int i=0;i<n-1;i++) {
            int j = urand48() % (n-i) + i;
            if(i!=j) swap(narray[i],narray[j]);
        }
    }
#endif

  float CInstanceP2D::get_param(int i) {
    switch (i) {
      case 0:
        return translation[0];
      case 1:
        return translation[1];
      case 2:
        return angle;
      case 3:
        return scale;
      default:
        throw "parameter index out of range";
    }
  }

  CInstanceP2D::CInstanceP2D() { init(); }

  void CInstanceP2D::init() {
    image_size = 512;
    model_size = 100;
    nclutter = igetenv("nclutter", 50);
    nmodel_total = 20;
    nmodel_unoccluded = 10;
    error = 5.0;
    aerror = 0.1;
  }

  
  /**
   * @author Shawn Le 
   * @date   07/Mar/2016
   * @brief  read input data instead of using generated random data
   */ 
  //int CInstanceP2D::readInputs(cv::Mat& img_out)
  int CInstanceP2D::readInputs()
  {
    int ret = 0;
    
    cv::Mat modelimg, imageimg; 
    modelimg = cv::imread("model.bmp",CV_LOAD_IMAGE_COLOR);
    if (modelimg.empty())
    {
      std::cout << "error! model.bmp cannot retrieve" << std::endl;
      ret = 1;
      return ret;
    }
    cv::resize(modelimg, modelimg, cv::Size(), 0.125, .125); 
    
    imageimg = cv::imread("image.bmp",CV_LOAD_IMAGE_COLOR);
    if (imageimg.empty())
    {
      std::cout << "error! image.bmp cannot retrieve" << std::endl;
      ret = 1;
      return ret;
    }
    cv::resize(imageimg, imageimg, cv::Size(), 0.125, .125); 
    imageimg.copyTo(disp_img);
    
    EdgeDetector ed; 
    EdgeDetector::EdgeData modelEdDat, imageEdDat;
    cv::Mat modelEdge, imageEdge;
    ed.EdgeDetect(modelimg, modelEdDat, modelEdge);
    ed.EdgeDetect(imageimg, imageEdDat, imageEdge);
    
    cv::imshow("model",modelEdge);
    cv::imshow("image",imageEdge);
    cv::waitKey(0);
 
    // pass inputs to compatible data format
    msources.clear();
    ipoints.clear();
    nmodel_total = modelEdDat.edgePtsOrg.size();
    
    printf("\n ### number of model's edge points = %d\n",nmodel_total);
    printf("\n ### number of image's edge points = %d\n",imageEdDat.edgePtsOrg.size());
    
    for (int i = 0; i < nmodel_total; i++) {
      Msource &m = msources.push();
      m.p = vec2((float)modelEdDat.edgePtsOrg[i].x, (float)modelEdDat.edgePtsOrg[i].y);
      //printf("M (%.2f, %.2f) ",m.p[0],m.p[1]);
      m.a = (float)atan2(modelEdDat.DyArr[i], modelEdDat.DxArr[i]);
    }
    for (int i = 0; i < imageEdDat.edgePtsOrg.size(); i++) {
      Ipoint &p = ipoints.push();
      p.p = vec2((float)imageEdDat.edgePtsOrg[i].x, (float)imageEdDat.edgePtsOrg[i].y);
      //printf("I (%.2f, %.2f) ",p.p[0],p.p[1]);
      p.a = (float)atan2(imageEdDat.DyArr[i], imageEdDat.DxArr[i]);
    }
        
    return ret;
  }
  
  
  /// generate a random test data i.e. model and image
  /**
   * @details 
   * \n 1) model points are randomly created i.e.
   * \n --> m.p vec2(urand(-model_size, model_size), urand(-model_size, model_size))
   * \n --> m.a = urand(0.0, 2 * M_PI)
   * \n 2) projected image points + noise are created 
   * \n --> p.p = cmul(rotation, msources[i].p) + translation + randomUniformVectorFromCircle(error);
   * \n --> p.a = msources[i].a + angle + urand(-aerror, aerror);
   * \n 3) cluttered image points are randomly created
   * @param[out] msources with msources[i].p and msources[i].a
   * @param[out] ipoints  with ipoints[i].p and ipoints[i].a
   */
  void CInstanceP2D::generate() {
    float dx = urand(0.0, image_size);
    float dy = urand(0.0, image_size);
    translation = vec2(dx, dy);
    angle = urand(0.0, 2.0 * M_PI);
    scale = urand(minscale, maxscale);
    vec2 rotation = vec2(scale * cos(angle), scale * sin(angle));
    msources.clear();
    ipoints.clear();
    for (int i = 0; i < nmodel_total; i++) {
      Msource &m = msources.push();
      m.p =
          vec2(urand(-model_size, model_size), urand(-model_size, model_size));
      m.a = urand(0.0, 2 * M_PI);
    }
    for (int i = 0; i < nmodel_unoccluded; i++) {
      Ipoint &p = ipoints.push();
      p.p = cmul(rotation, msources[i].p) + translation +
            randomUniformVectorFromCircle(error);
      p.a = msources[i].a + angle + urand(-aerror, aerror);
    }
    shuffle(msources);
    for (int i = 0; i < nclutter; i++) {
      Ipoint &p = ipoints.push();
      p.p = vec2(urand(0, image_size), urand(0, image_size));
      p.a = urand(0.0, 2 * M_PI);
    }
    shuffle(ipoints);
  }

  void CInstanceP2D::set_image_size(int r) { image_size = r; }
  void CInstanceP2D::set_model_size(int r) { model_size = r; }
  void CInstanceP2D::set_nclutter(int v) { nclutter = v; }
  void CInstanceP2D::set_nmodel_total(int v) { nmodel_total = v; }
  void CInstanceP2D::set_nmodel_unoccluded(int v) { nmodel_unoccluded = v; }
  void CInstanceP2D::set_error(float v) { error = v; }
  void CInstanceP2D::set_aerror(float v) { aerror = v; }
  void CInstanceP2D::set_srange(float min, float max) {
    minscale = min;
    maxscale = max;
  }
  int CInstanceP2D::nimage() { return ipoints.length(); }
  /// @param[out] a -> orientation value e.g. angle of the point
  void CInstanceP2D::get_image(float &x, float &y, float &a, int i) {
    x = ipoints[i].p[0];
    y = ipoints[i].p[1];
    a = ipoints[i].a;
  }
  int CInstanceP2D::nmodel() { return msources.length(); }
  void CInstanceP2D::get_model(float &x, float &y, float &a, int i) {
    x = msources[i].p[0];
    y = msources[i].p[1];
    a = msources[i].a;
  }

  CInstanceP2D::~CInstanceP2D() {}
  
}


//InstanceP2D *makeInstanceP2D() { return new lumo_cinstancep2d::InstanceP2D(); }
/**
 * @author  Shawn Le
 * @details replace his prototype returning a parent class to returning a child class, allowing new implementations of child class can be used.
 *          This also requires child class's prototype should be moved to header for global calls and uses
 */
lumo_cinstancep2d::CInstanceP2D *makeInstanceP2D() { return new lumo_cinstancep2d::CInstanceP2D(); }

/**
 * @author  Shawn Le
 * @brief   I create another prototype based on makeInstanceP2D so that new makeCInstanceP2D can return my new class CInstanceP2D  
 * @details his makeInstanceP2D calls implementation of CInstanceP2D but returns a InstanceP2D which is the parent class and is a virtual class
 *          it makes new implementations of CInstanceP2D is not acceptable by compiler since they don't exist in the parent class
 */
lumo_cinstancep2d::CInstanceP2D *makeCInstanceP2D() { return new lumo_cinstancep2d::CInstanceP2D(); }
