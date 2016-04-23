/**
 * @file   EdgeDetector.cpp 
 * @brief  EdgeDetector in-sight code 
 * @author Son Le (Shawn)
 * @date   15/10/2015
 */

#include "EdgeDetector.hpp"
#include <boost/graph/graph_concepts.hpp>

const int DEBUG_MODE = 0;

using namespace std;
using namespace cv;


EdgeDetector::EdgeDetector()
{
}




EdgeDetector::~EdgeDetector()
{
}




int EdgeDetector::EdgeDetect(const Mat image, EdgeData& edDat, Mat& edge_ou)
{
  int ret = 0;
  
  image.copyTo(img);
  
  // **** to gray
  if (img.type() == CV_8UC3)
  {
    cvtColor(img, imgGray, CV_BGR2GRAY);   
  } else
  {
    cout << "[EdgeDetector][ERROR] image must be colored" << endl; 
    return 1;
  }
  
  
  // **** filter
  //GaussianBlur(imgGray, imgGray, Size(3,3), 0);
  
  
  // **** gradient
  // Sobel X -> Sobel(src, dst, ddepth, dx, dy, ksize, scale, delta, border_type)
  // note: where dx is column-wise, dy is row-wise
  Sobel(imgGray, Sx, CV_16SC1, 1, 0, 3); 
  Sobel(imgGray, Sy, CV_16SC1, 0, 1, 3); 
  
  if (DEBUG_MODE)
  {
    imshow("sb x", Sx);
    waitKey(0);
    // Sobel Y
    imshow("sb y", Sy);
    waitKey(0);
  }
  
  
  // **** Phase and magnitude calculation
    const short* _sdx; 
    const short* _sdy; 
    double fdx,fdy; 
    double MagG, DirG;
    double MaxGradient=-99999.99;
    double direction;
    vector<int> orients;
    int count = 0,i,j; // count variable;
    
    Mat magMat = Mat::zeros(Sx.rows, Sx.cols, CV_64F);
    
    for( i = 1; i < Sx.rows-1; i++ )
    {
        for( j = 1; j < Sx.cols-1; j++ )
        {        
                fdx = (float) Sx.at<short>(i,j); 
                fdy = (float) Sy.at<short>(i,j);

                MagG = sqrt((float)(fdx*fdx) + (float)(fdy*fdy)); //Magnitude = Sqrt(gx^2 +gy^2)
                direction =cvFastArctan((float)fdy,(float)fdx);  //Direction = invtan (Gy / Gx)
                magMat.at<double>(i,j) = MagG;
                
                if(MagG>MaxGradient)
                    MaxGradient=MagG; // get maximum gradient value for normalizing.

                    // get closest angle from 0, 45, 90, 135 set
                        if ( (direction>0 && direction < 22.5) || (direction >157.5 && direction < 202.5) || (direction>337.5 && direction<360)  )
                            direction = 0;
                        else if ( (direction>22.5 && direction < 67.5) || (direction >202.5 && direction <247.5)  )
                            direction = 45;
                        else if ( (direction >67.5 && direction < 112.5)||(direction>247.5 && direction<292.5) )
                            direction = 90;
                        else if ( (direction >112.5 && direction < 157.5)||(direction>292.5 && direction<337.5) )
                            direction = 135;
                        else 
                            direction = 0;
                
            orients.push_back((int)direction);
            count++;
        }
    }
  
  
  count=0; // init count
  
  // **** non maximum suppression
  Mat nms;
  nms = Mat::zeros(Sx.rows, Sx.cols, CV_32F); 
  double leftPixel,rightPixel;
    
    for( i = 1; i < Sx.rows-1; i++ )
    {
        for( j = 1; j < Sx.cols-1; j++ )
        {
                switch ( orients[count] )
                {
                   case 0:
                        leftPixel  = magMat.at<double>(i,j-1); 
                        rightPixel = magMat.at<double>(i,j+1);
                        break;
                    case 45:
                        leftPixel  = magMat.at<double>(i-1,j+1);
                        rightPixel = magMat.at<double>(i+1,j-1);
                        break;
                    case 90:
                        leftPixel  = magMat.at<double>(i-1,j);
                        rightPixel = magMat.at<double>(i+1,j);
                        break;
                    case 135:
                        leftPixel  = magMat.at<double>(i-1,j-1);
                        rightPixel = magMat.at<double>(i+1,j+1);
                        break;
                 }
                // compare current pixels value with adjacent pixels
                if (( magMat.at<double>(i,j) < leftPixel ) || (magMat.at<double>(i,j) < rightPixel ) )
                  nms.at<float>(i,j) = 0;
                else
                  nms.at<float>(i,j) = (uchar)(magMat.at<double>(i,j)/MaxGradient*255.);
            
                count++;
            }
        }
  
  
  Mat nmsMat(nms), nmsMatdsp; 
  nmsMat.copyTo(nmsMatdsp);
  threshold(nmsMatdsp, nmsMatdsp, 0, 255, CV_THRESH_BINARY);
  
  if (DEBUG_MODE)
  {
    cout << nmsMatdsp.rows << endl;
    cout << nmsMatdsp.cols << endl;
    imshow("nms edge", nmsMatdsp);
    waitKey(0);
  }
  
  
  int RSum=0,CSum=0;
  int curX,curY;
  int flag=1;

  // **** Hysterisis threshold
  double maxContrast = 10.;
  double minContrast = 100.;
  vector<Point2i> coords;
  vector<double> edgeDx, edgeDy, edgeMag;
  int noOfCoords = 0; 
  
    for( i = 1; i < Sx.rows-1; i++ )
    {
        for( j = 1; j < Sx.cols; j++ )
        {
            fdx = (float) Sx.at<short>(i,j); 
            fdy = (float) Sy.at<short>(i,j);
                
            MagG = sqrt(fdx*fdx + fdy*fdy); //Magnitude = Sqrt(gx^2 +gy^2)
            DirG =cvFastArctan((float)fdy,(float)fdx);   //Direction = tan(y/x)
        
            flag=1;
            if((double)(nms.at<float>(i,j)) < maxContrast)
            {
                if((double)(nms.at<float>(i,j))< minContrast)
                {
                    nms.at<float>(i,j) = 0;
                    flag=0; // remove from edge
                }
                else
                {   // if any of 8 neighboring pixel is not greater than max contraxt remove from edge
                    if( ((double)(nms.at<float>(i-1,j-1) < maxContrast))    &&
                        ((double)(nms.at<float>(i-1,j)   < maxContrast))    &&
                        ((double)(nms.at<float>(i-1,j+1) < maxContrast))    &&
                        ((double)(nms.at<float>(i,j-1)   < maxContrast))    &&
                        ((double)(nms.at<float>(i,j+1)   < maxContrast))    &&
                        ((double)(nms.at<float>(i+1,j-1) < maxContrast))    &&
                        ((double)(nms.at<float>(i+1,j)   < maxContrast))    &&
                        ((double)(nms.at<float>(i+1,j+1) < maxContrast))       )
                    {
                        nms.at<float>(i,j) = 0;
                        flag=0;
                    }
                }
                
            }
            
            // save selected edge information
            // curX=i; curY=j;  // for consistency, X-> cols Y-> rows
            curX=j; curY=i;
            if(flag!=0)
            {
                if(fdx!=0 || fdy!=0)
                {       
                    RSum=RSum+curX; CSum=CSum+curY; // Row sum and column sum for center of gravity
                    
                    Point2i p; 
                    p.x = curX; p.y = curY;
                    coords.push_back(p);
                    edgeDx.push_back(fdx);
                    edgeDy.push_back(fdy);
                    
                    //handle divide by zero
                    if(MagG!=0)
                        edgeMag.push_back(1/MagG);  // gradient magnitude 
                    else
                        edgeMag.push_back(0);
                                                            
                    noOfCoords++;
                }
            }
        }
    }

    Point2i COG;
    COG.x = RSum/noOfCoords; // center of gravity
    COG.y = CSum/noOfCoords; // center of gravity
        
    // change coordinates to reflect center of gravity
    edDat.edgePtsOrg = coords; 
    for(int m=0;m<noOfCoords ;m++)
    {
        int temp;

        temp = coords[m].x;
        coords[m].x = temp-COG.x;
        temp = coords[m].y;
        coords[m].y = temp-COG.y;
    }
  
  
  
  // **** edge linking
  // RECOMMENDATION: 0.8 <= tHigh <= 0.9,  0.3 <= tLow <= 0.5
  
  Sx.copyTo(edDat.DxMat);
  Sy.copyTo(edDat.DyMat);
  magMat.copyTo(edDat.magMat);
  
  edDat.DxArr = edgeDx;
  edDat.DyArr = edgeDy;
  edDat.magArr =  edgeMag;
  edDat.edgePts = coords;  
  
  nms.copyTo(edge_ou);
  
  return ret;
}
