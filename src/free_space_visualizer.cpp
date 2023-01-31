//
// Created by Jacobus Conradi on 04.02.22.
//

#include "free_space_visualizer.h"


#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

using std::cout;
using std::endl;
using namespace cv;

FreeSpaceVisualizer::FreeSpaceVisualizer(FreeSpace &freespace) : freespace(freespace) {
}
//circleradius
#define CR 10

void FreeSpaceVisualizer::show(bool withPoints){

    if(freespace.xSize() == 0)
        return;

    int CS = 8;



    int scrolHight = 0;
    int scrolWidth = 0;


    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img((freespace.ySize() + 2) * CS, (freespace.xSize() + 2) * CS, CV_8UC3, Scalar(255, 255, 255));

        //draw vertexlines
        for (int i = 0; i <= freespace.xSize(); ++i) {
            int x = (i + 1) * CS;
            int y1 = CS;
            int y2 = (freespace.ySize() + 1) * CS;
            line(img, Point2i(x, y1), Point2i(x, y2), Scalar(0, 0, 0), 3);
        }

        for (int i = 0; i <= freespace.ySize(); ++i) {
            int y = (i + 1) * CS;
            int x1 = CS;
            int x2 = (freespace.xSize() + 1) * CS;
            line(img, Point2i(x1, y), Point2i(x2, y), Scalar(0, 0, 0), 3);
        }

        for (int x = 0; x < freespace.xSize(); ++x) {
            for (int y = 0; y < freespace.ySize(); ++y) {

                Cell *cell = freespace.getCell(x, y);

                if (cell->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //eight points, eight lines

                Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 3);
                line(img, p2, p3, Scalar(0, 0, 255), 3);
                line(img, p3, p4, Scalar(0, 0, 255), 3);
                line(img, p4, p5, Scalar(0, 0, 255), 3);
                line(img, p5, p6, Scalar(0, 0, 255), 3);
                line(img, p6, p7, Scalar(0, 0, 255), 3);
                line(img, p7, p8, Scalar(0, 0, 255), 3);
                line(img, p8, p1, Scalar(0, 0, 255), 3);

/*
                for (auto p: cell->importantUpYs) {
                    Point2d ip = Point2d(p.first.x, -p.first.y) * CS + origin;
                    circle(img, ip, CR, Scalar(0, 255, 0), -1);
                }
                for (auto p: cell->importantDownYs) {
                    Point2d ip = Point2d(p.first.x, -p.first.y) * CS + origin;
                    circle(img, ip, CR, Scalar(255, 0, 0), -1);
                }
                */
            }
        }

        //draw dots
        if(withPoints){
            for (auto p : freespace.pointvector){

                auto x = p.x;
                auto y = p.y;
                Point2d originA = Point2d((x.id + 1) * CS, (freespace.ySize() + 1 - y.id) * CS);
                Point2d point = Point2d(x.t, -y.t) * CS + originA;
                circle(img,point, 5, Scalar(0, 255, 0), -1);
            }
        }

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,512);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,4);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

void FreeSpaceVisualizer::show(FreeSpacePoint &a, FreeSpacePoint &b) {

    if(freespace.xSize() == 0)
        return;

    int CS = 64;



    int scrolHight = 0;
    int scrolWidth = 0;


    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img((freespace.ySize() + 2) * CS, (freespace.xSize() + 2) * CS, CV_8UC3, Scalar(255, 255, 255));

        //draw vertexlines
        for (int i = 0; i <= freespace.xSize(); ++i) {
            int x = (i + 1) * CS;
            int y1 = CS;
            int y2 = (freespace.ySize() + 1) * CS;
            line(img, Point2i(x, y1), Point2i(x, y2), Scalar(0, 0, 0), 3);
        }

        for (int i = 0; i <= freespace.ySize(); ++i) {
            int y = (i + 1) * CS;
            int x1 = CS;
            int x2 = (freespace.xSize() + 1) * CS;
            line(img, Point2i(x1, y), Point2i(x2, y), Scalar(0, 0, 0), 3);
        }

        for (int x = 0; x < freespace.xSize(); ++x) {
            for (int y = 0; y < freespace.ySize(); ++y) {

                Cell *cell = freespace.getCell(x, y);

                if (cell->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //eight points, eight lines

                Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 3);
                line(img, p2, p3, Scalar(0, 0, 255), 3);
                line(img, p3, p4, Scalar(0, 0, 255), 3);
                line(img, p4, p5, Scalar(0, 0, 255), 3);
                line(img, p5, p6, Scalar(0, 0, 255), 3);
                line(img, p6, p7, Scalar(0, 0, 255), 3);
                line(img, p7, p8, Scalar(0, 0, 255), 3);
                line(img, p8, p1, Scalar(0, 0, 255), 3);
            }
        }

        //A

        int imageYA = (freespace.ySize() + 1 - a.y.id) * CS;
        Point2d originA = Point2d((a.x.id + 1) * CS, imageYA);
        Point2d p1 = Point2d(a.x.t, -a.y.t) * CS + originA;

        int imageYB = (freespace.ySize() + 1 - b.y.id) * CS;
        Point2d originB = Point2d((b.x.id + 1) * CS, imageYB);
        Point2d p2 = Point2d(b.x.t, -b.y.t) * CS + originB;

        line(img,p1 , p2, Scalar(0,255,0),5);

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,512);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,16);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}


void FreeSpaceVisualizer::show(FreeSpacePoint &a, FreeSpacePoint &b, ParamPoint x1, ParamPoint x2) {

    if(freespace.xSize() == 0)
        return;

    int CS = 64;



    int scrolHight = 0;
    int scrolWidth = 0;


    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img((freespace.ySize() + 2) * CS, (freespace.xSize() + 2) * CS, CV_8UC3, Scalar(255, 255, 255));

        //draw vertexlines
        for (int i = 0; i <= freespace.xSize(); ++i) {
            int x = (i + 1) * CS;
            int y1 = CS;
            int y2 = (freespace.ySize() + 1) * CS;
            line(img, Point2i(x, y1), Point2i(x, y2), Scalar(0, 0, 0), 3);
        }

        for (int i = 0; i <= freespace.ySize(); ++i) {
            int y = (i + 1) * CS;
            int x1 = CS;
            int x2 = (freespace.xSize() + 1) * CS;
            line(img, Point2i(x1, y), Point2i(x2, y), Scalar(0, 0, 0), 3);
        }

        for (int x = 0; x < freespace.xSize(); ++x) {
            for (int y = 0; y < freespace.ySize(); ++y) {

                Cell *cell = freespace.getCell(x, y);

                if (cell->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //eight points, eight lines

                Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 3);
                line(img, p2, p3, Scalar(0, 0, 255), 3);
                line(img, p3, p4, Scalar(0, 0, 255), 3);
                line(img, p4, p5, Scalar(0, 0, 255), 3);
                line(img, p5, p6, Scalar(0, 0, 255), 3);
                line(img, p6, p7, Scalar(0, 0, 255), 3);
                line(img, p7, p8, Scalar(0, 0, 255), 3);
                line(img, p8, p1, Scalar(0, 0, 255), 3);
            }
        }


        //startline

        //X1
        int x = (freespace.xSize() + 2) * CS;
        int imageYX1 = (freespace.ySize() + 1 - x1.id) * CS;
        Point2d originX11 = Point2d(0, imageYX1);
        Point2d originX12 = Point2d(x, imageYX1);
        Point2d x11 = Point2d(0, -x1.t) * CS + originX11;
        Point2d x12 = Point2d(0, -x1.t) * CS + originX12;

        line(img,x11 , x12, Scalar(255,0,0),2);


        //X2
        int imageYX2 = (freespace.ySize() + 1 - x2.id) * CS;
        Point2d originX21 = Point2d(0, imageYX2);
        Point2d originX22 = Point2d(x, imageYX2);
        Point2d x21 = Point2d(0, -x2.t) * CS + originX21;
        Point2d x22 = Point2d(0, -x2.t) * CS + originX22;

        line(img,x21 , x22, Scalar(255,0,0),2);

        //A
        int imageYA = (freespace.ySize() + 1 - a.y.id) * CS;
        Point2d originA = Point2d((a.x.id + 1) * CS, imageYA);
        Point2d p1 = Point2d(a.x.t, -a.y.t) * CS + originA;

        //B
        int imageYB = (freespace.ySize() + 1 - b.y.id) * CS;
        Point2d originB = Point2d((b.x.id + 1) * CS, imageYB);
        Point2d p2 = Point2d(b.x.t, -b.y.t) * CS + originB;

        line(img,p1 , p2, Scalar(0,255,0),5);

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,512);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,16);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

void FreeSpaceVisualizer::showCandidate(Candidate & c, int curveId){
    if(freespace.xSize() == 0)
        return;

    int CS = 64;



    int scrolHight = 0;
    int scrolWidth = 0;


    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img((freespace.ySize() + 2) * CS, (freespace.xSize() + 2) * CS, CV_8UC3, Scalar(255, 255, 255));

        //draw vertexlines
        for (int i = 0; i <= freespace.xSize(); ++i) {
            int x = (i + 1) * CS;
            int y1 = CS;
            int y2 = (freespace.ySize() + 1) * CS;
            line(img, Point2i(x, y1), Point2i(x, y2), Scalar(0, 0, 0), 3);
        }

        for (int i = 0; i <= freespace.ySize(); ++i) {
            int y = (i + 1) * CS;
            int x1 = CS;
            int x2 = (freespace.xSize() + 1) * CS;
            line(img, Point2i(x1, y), Point2i(x2, y), Scalar(0, 0, 0), 3);
        }

        for (int x = 0; x < freespace.xSize(); ++x) {
            for (int y = 0; y < freespace.ySize(); ++y) {

                Cell *cell = freespace.getCell(x, y);

                if (cell->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //eight points, eight lines

                Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 3);
                line(img, p2, p3, Scalar(0, 0, 255), 3);
                line(img, p3, p4, Scalar(0, 0, 255), 3);
                line(img, p4, p5, Scalar(0, 0, 255), 3);
                line(img, p5, p6, Scalar(0, 0, 255), 3);
                line(img, p6, p7, Scalar(0, 0, 255), 3);
                line(img, p7, p8, Scalar(0, 0, 255), 3);
                line(img, p8, p1, Scalar(0, 0, 255), 3);
            }
        }


        //startline
        ParamPoint yStart = c.getStart();
        ParamPoint yEnd = c.getEnd();

        //X1
        int x = (freespace.xSize() + 2) * CS;
        int imageYX1 = (freespace.ySize() + 1 - yStart.id) * CS;
        Point2d originX11 = Point2d(0, imageYX1);
        Point2d originX12 = Point2d(x, imageYX1);
        Point2d x11 = Point2d(0, -yStart.t) * CS + originX11;
        Point2d x12 = Point2d(0, -yStart.t) * CS + originX12;

        line(img, x11, x12, Scalar(255, 0, 0), 2);


        //X2
        int imageYX2 = (freespace.ySize() + 1 - yEnd.id) * CS;
        Point2d originX21 = Point2d(0, imageYX2);
        Point2d originX22 = Point2d(x, imageYX2);
        Point2d x21 = Point2d(0, -yEnd.t) * CS + originX21;
        Point2d x22 = Point2d(0, -yEnd.t) * CS + originX22;

        line(img, x21, x22, Scalar(255, 0, 0), 2);



        for(auto cov : c.visualMatchings) {
            if(cov.curveIdx != curveId)
                continue;
            ParamPoint xStart = cov.getStart();
            ParamPoint xEnd = cov.getEnd();

            FreeSpacePoint a = {xStart,yStart};
            FreeSpacePoint b = {xEnd,yEnd};
            //A
            int imageYA = (freespace.ySize() + 1 - a.y.id) * CS;
            Point2d originA = Point2d((a.x.id + 1) * CS, imageYA);
            Point2d p1 = Point2d(a.x.t, -a.y.t) * CS + originA;

            //B
            int imageYB = (freespace.ySize() + 1 - b.y.id) * CS;
            Point2d originB = Point2d((b.x.id + 1) * CS, imageYB);
            Point2d p2 = Point2d(b.x.t, -b.y.t) * CS + originB;

            line(img, p1, p2, Scalar(0, 255, 0), 5);
        }

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,512);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,16);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

void FreeSpaceVisualizer::showCandidates(std::vector<Candidate> candidates, int curveId){
    if(freespace.xSize() == 0)
        return;

    int CS = 8;



    int scrolHight = 0;
    int scrolWidth = 0;


    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img((freespace.ySize() + 2) * CS, (freespace.xSize() + 2) * CS, CV_8UC3, Scalar(255, 255, 255));

        //draw vertexlines
        for (int i = 0; i <= freespace.xSize(); ++i) {
            int x = (i + 1) * CS;
            int y1 = CS;
            int y2 = (freespace.ySize() + 1) * CS;
            line(img, Point2i(x, y1), Point2i(x, y2), Scalar(0, 0, 0), 3);
        }

        for (int i = 0; i <= freespace.ySize(); ++i) {
            int y = (i + 1) * CS;
            int x1 = CS;
            int x2 = (freespace.xSize() + 1) * CS;
            line(img, Point2i(x1, y), Point2i(x2, y), Scalar(0, 0, 0), 3);
        }

        for (int x = 0; x < freespace.xSize(); ++x) {
            for (int y = 0; y < freespace.ySize(); ++y) {

                Cell *cell = freespace.getCell(x, y);

                if (cell->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //eight points, eight lines

                Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 3);
                line(img, p2, p3, Scalar(0, 0, 255), 3);
                line(img, p3, p4, Scalar(0, 0, 255), 3);
                line(img, p4, p5, Scalar(0, 0, 255), 3);
                line(img, p5, p6, Scalar(0, 0, 255), 3);
                line(img, p6, p7, Scalar(0, 0, 255), 3);
                line(img, p7, p8, Scalar(0, 0, 255), 3);
                line(img, p8, p1, Scalar(0, 0, 255), 3);
            }
        }





        for(auto c : candidates) {

            //startline
            ParamPoint yStart = c.getStart();
            ParamPoint yEnd = c.getEnd();

            //X1
            int x = (freespace.xSize() + 2) * CS;
            int imageYX1 = (freespace.ySize() + 1 - yStart.id) * CS;
            Point2d originX11 = Point2d(0, imageYX1);
            Point2d originX12 = Point2d(x, imageYX1);
            Point2d x11 = Point2d(0, -yStart.t) * CS + originX11;
            Point2d x12 = Point2d(0, -yStart.t) * CS + originX12;

            line(img, x11, x12, Scalar(255, 0, 0), 2);


            //X2
            int imageYX2 = (freespace.ySize() + 1 - yEnd.id) * CS;
            Point2d originX21 = Point2d(0, imageYX2);
            Point2d originX22 = Point2d(x, imageYX2);
            Point2d x21 = Point2d(0, -yEnd.t) * CS + originX21;
            Point2d x22 = Point2d(0, -yEnd.t) * CS + originX22;

            line(img, x21, x22, Scalar(255, 0, 0), 2);

            for (auto cov: c.visualMatchings) {
                if (cov.curveIdx != curveId)
                    continue;
                ParamPoint xStart = cov.getStart();
                ParamPoint xEnd = cov.getEnd();

                FreeSpacePoint a = {xStart, yStart};
                FreeSpacePoint b = {xEnd, yEnd};
                //A
                int imageYA = (freespace.ySize() + 1 - a.y.id) * CS;
                Point2d originA = Point2d((a.x.id + 1) * CS, imageYA);
                Point2d p1 = Point2d(a.x.t, -a.y.t) * CS + originA;

                //B
                int imageYB = (freespace.ySize() + 1 - b.y.id) * CS;
                Point2d originB = Point2d((b.x.id + 1) * CS, imageYB);
                Point2d p2 = Point2d(b.x.t, -b.y.t) * CS + originB;

                line(img, p1, p2, Scalar(0, 255, 0), 5);
            }
        }

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,512);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,16);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

int getRealY(int y, int maxY){
    return (maxY - y);
}


void FreeSpacesVisualizer::showCandidates(CandidateSet &candidates, std::vector<std::pair<int, int>>& indices) {
    //this will be awful
    int CS = 4;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<int> partialYSums;
    int ySize = 1;
    for(auto col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
    }

    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto fs : freespaces[0]){
        partialXSums.push_back(xSize);
        xSize += (fs.xSize() + 2);
    }

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(xSize*CS, ySize*CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                //freespaceorigin
                int xOrigin = partialXSums[xi];
                int yOrigin = partialYSums[yi];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code
                FreeSpace &freespace = freespaces[yi][xi];

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1) + fsOrigin, Point2d(x, -y2) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int x = 0; x < freespace.xSize(); ++x) {
                    for (int y = 0; y < freespace.ySize(); ++y) {

                        Cell *cell = freespace.getCell(x, y);

                        if (cell->is_empty())
                            continue;

                        //stupid mapping only flip y
                        int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //eight points, eight lines

                        Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 1);
                        line(img, p2, p3, Scalar(0, 0, 255), 1);
                        line(img, p3, p4, Scalar(0, 0, 255), 1);
                        line(img, p4, p5, Scalar(0, 0, 255), 1);
                        line(img, p5, p6, Scalar(0, 0, 255), 1);
                        line(img, p6, p7, Scalar(0, 0, 255), 1);
                        line(img, p7, p8, Scalar(0, 0, 255), 1);
                        line(img, p8, p1, Scalar(0, 0, 255), 1);
                    }
                }
            }
        }

        //candidates

        for(auto i:indices){
            int baseCurveIdx = i.first;
            //first draw horizontal lines
            Candidate candidate = candidates.getCandidate(i);
            double y1 = ySize*CS - (candidate.getStart().t + candidate.getStart().id + partialYSums[baseCurveIdx]+1)*CS;
            double y2 = ySize*CS - (candidate.getEnd().t + candidate.getEnd().id + partialYSums[baseCurveIdx]+1)*CS;
            line(img, Point2d(0,y1), Point2d(xSize*CS,y1), Scalar(255, 0, 0), 1);
            line(img, Point2d(0,y2), Point2d(xSize*CS,y2), Scalar(255, 0, 0), 1);

            for(auto cov : candidate.visualMatchings){
                int targetCurveIdx = cov.curveIdx;
                double x1 = (cov.getStart().t + cov.getStart().id + partialXSums[targetCurveIdx]+1)*CS;
                double x2 = (cov.getEnd().t + cov.getEnd().id + partialXSums[targetCurveIdx]+1)*CS;

                line(img, Point2d(x1,y1), Point2d(x2,y2), Scalar(0, 255, 0), 3);
                
            }


        }


        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,128);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,1);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

void FreeSpacesVisualizer::showCandidates(std::vector<Candidate> candidates) {
    //this will be awful
    int CS = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<int> partialYSums;
    int ySize = 1;
    for(auto col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
    }

    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto fs : freespaces[0]){
        partialXSums.push_back(xSize);
        xSize += (fs.xSize() + 2);
    }

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(xSize*CS, ySize*CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                //freespaceorigin
                int xOrigin = partialXSums[xi];
                int yOrigin = partialYSums[yi];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code
                FreeSpace &freespace = freespaces[yi][xi];

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1) + fsOrigin, Point2d(x, -y2) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int x = 0; x < freespace.xSize(); ++x) {
                    for (int y = 0; y < freespace.ySize(); ++y) {

                        Cell *cell = freespace.getCell(x, y);

                        if (cell->is_empty())
                            continue;

                        //stupid mapping only flip y
                        int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //eight points, eight lines

                        Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 1);
                        line(img, p2, p3, Scalar(0, 0, 255), 1);
                        line(img, p3, p4, Scalar(0, 0, 255), 1);
                        line(img, p4, p5, Scalar(0, 0, 255), 1);
                        line(img, p5, p6, Scalar(0, 0, 255), 1);
                        line(img, p6, p7, Scalar(0, 0, 255), 1);
                        line(img, p7, p8, Scalar(0, 0, 255), 1);
                        line(img, p8, p1, Scalar(0, 0, 255), 1);
                    }
                }
            }
        }

        //candidates

        for(int i=0;i<candidates.size();++i){
            auto candidatepair = candidates[i];
            bool last = i == candidates.size()-1;
            int baseCurveIdx = candidatepair.getIndex();
            //first draw horizontal lines
            Candidate candidate = candidatepair;
            double y1 = ySize*CS - (candidate.getStart().t + candidate.getStart().id + partialYSums[baseCurveIdx]+1)*CS;
            double y2 = ySize*CS - (candidate.getEnd().t + candidate.getEnd().id + partialYSums[baseCurveIdx]+1)*CS;
            line(img, Point2d(0,y1), Point2d(xSize*CS,y1), Scalar(255, 0, 0), 1);
            line(img, Point2d(0,y2), Point2d(xSize*CS,y2), Scalar(255, 0, 0), 1);

            for(auto cov : candidate.visualMatchings){
                int targetCurveIdx = cov.getIndex();
                double x1 = (cov.getStart().t + cov.getStart().id + partialXSums[targetCurveIdx]+1)*CS;
                double x2 = (cov.getEnd().t + cov.getEnd().id + partialXSums[targetCurveIdx]+1)*CS;

                line(img, Point2d(x1,y1), Point2d(x2,y2), Scalar(last?255:0, last?0:255, 0), 3);

            }


        }


        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,128);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,1);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

void FreeSpacesVisualizer::show() {
//this will be awful
    int CS = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<int> partialYSums;
    int ySize = 1;
    for(auto col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
    }

    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto fs : freespaces[0]){
        partialXSums.push_back(xSize);
        xSize += (fs.xSize() + 2);
    }

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(xSize*CS, ySize*CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                //freespaceorigin
                int xOrigin = partialXSums[xi];
                int yOrigin = partialYSums[yi];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code
                FreeSpace &freespace = freespaces[yi][xi];

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1) + fsOrigin, Point2d(x, -y2) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y) + fsOrigin, Scalar(0, 0, 0), 1);
                }

                for (int x = 0; x < freespace.xSize(); ++x) {
                    for (int y = 0; y < freespace.ySize(); ++y) {

                        Cell *cell = freespace.getCell(x, y);

                        if (cell->is_empty())
                            continue;

                        //stupid mapping only flip y
                        int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //eight points, eight lines

                        Point2d p1 = Point2d(cell->leftPair.first.x, -cell->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cell->leftPair.second.x, -cell->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cell->topPair.first.x, -cell->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cell->topPair.second.x, -cell->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cell->rightPair.first.x, -cell->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cell->rightPair.second.x, -cell->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cell->bottomPair.first.x, -cell->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cell->bottomPair.second.x, -cell->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 1);
                        line(img, p2, p3, Scalar(0, 0, 255), 1);
                        line(img, p3, p4, Scalar(0, 0, 255), 1);
                        line(img, p4, p5, Scalar(0, 0, 255), 1);
                        line(img, p5, p6, Scalar(0, 0, 255), 1);
                        line(img, p6, p7, Scalar(0, 0, 255), 1);
                        line(img, p7, p8, Scalar(0, 0, 255), 1);
                        line(img, p8, p1, Scalar(0, 0, 255), 1);
                    }
                }
            }
        }

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,128);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,1);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}

Scalar _labelColor(Label l){
    switch (l) {

        case transition:
            return Scalar(0,0,0);
        case walk:
            return Scalar(0,255,255);
        case jump:
            return Scalar(64,190,255);
        case punch:
            return Scalar(128,128,255);
        case leg_kick:
            return Scalar(190,64,255);
        case squat:
            return Scalar(255,0,255);
        case run:
            return Scalar(255,64,190);
        case stand:
            return Scalar(255,128,128);
        case arm_up:
            return Scalar(255,190,64);
        case drink:
            return Scalar(255,255,0);
        case stretch:
            return Scalar(190,255,64);
        case slap:
            return Scalar(128,255,128);
        case turn:
            return Scalar(64,255,190);
        default:
            return Scalar(0,0,0);
    }
}

void ClusteringVisulaizer::showClustering(Curves c, std::vector<std::vector<std::pair<Label, ParamPoint>>> groundthruth,
                                          std::vector<Candidate> candidates) {

    int scale = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<Label> labeling(candidates.size());
    int labelingIdx = 0;

    std::vector<int> partialLengths;
    partialLengths.push_back(1);
    for (auto curve: c) {
        partialLengths.push_back(partialLengths[partialLengths.size() - 1] + 1 + curve.size());
    }
    int xSize = (partialLengths[partialLengths.size() - 1]);
    int ySize = 1 + (_labelend - _labelstart - 1) + 5 + candidates.size();

    std::cout << xSize << " x " << ySize << std::endl;

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(scale * ySize, scale * xSize, CV_8UC3, Scalar(255, 255, 255));

        for (int i = 0; i < c.size(); i++) {
            line(img, Point2d(partialLengths[i], 1) * scale, Point2d(partialLengths[i + 1]-1, 1) * scale, Scalar(0, 0, 0),
                 5);
            circle(img, Point2d(partialLengths[i], 1) * scale, 5, Scalar(0, 0, 0), -1);
            circle(img, Point2d(partialLengths[i + 1]-1, 1) * scale, 5, Scalar(0, 0, 0), -1);
        }

        for(Label label = (Label)(_labelstart+1);label!=_labelend;label = (Label)(label+1)){
            double y = 1+label;
            line(img, Point2d(partialLengths[0], y) * scale, Point2d(partialLengths[partialLengths.size()-1], y) * scale, Scalar(0, 0, 0),
                 1);
        }

        for (int curveindex = 0; curveindex < c.size(); curveindex++) {
            ParamPoint start = {0, 0};
            for (auto assignment: groundthruth[curveindex]) {
                double y = 1 + assignment.first;

                double intxs = partialLengths[curveindex] + start.id + start.t;
                double intxt = partialLengths[curveindex] + assignment.second.id + assignment.second.t;

                line(img, Point2d(intxs, y) * scale, Point2d(intxt, y) * scale, _labelColor(assignment.first), 5);
                circle(img, Point2d(intxs, y) * scale, 5, _labelColor(assignment.first), -1);
                circle(img, Point2d(intxt, y) * scale, 5, _labelColor(assignment.first), -1);

                start = assignment.second;
            }
        }

        for (int i = 0; i < candidates.size(); i++) {
            auto wrappercandidate = candidates[i];
            Candidate candidate = wrappercandidate;
            double y = 1 + (_labelend - _labelstart - 1) + 5 + i;
            line(img, Point2d(partialLengths[0], y) * scale, Point2d(partialLengths[partialLengths.size()-1], y) * scale, Scalar(0, 0, 0),
                 1);
            for (auto matching: candidate.visualMatchings) {
                double intxs = partialLengths[matching.curveIdx] + matching.start.id + matching.start.t;
                double intxt = partialLengths[matching.curveIdx] + matching.end.id + matching.end.t;

                line(img, Point2d(intxs, y+0.1) * scale, Point2d(intxt, y-0.1) * scale, _labelColor(labeling[i]), 5);
                circle(img, Point2d(intxs, y+0.1) * scale, 5, _labelColor(labeling[i]), -1);
                circle(img, Point2d(intxt, y-0.1) * scale, 5, _labelColor(labeling[i]), -1);
            }
        }
        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH)*scrolHight/1000;
            int truescrolWidth = (img.cols - winW)*scrolWidth/1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                scale = std::min(scale*2,128);
                redraw = true;
                break;
            }
            if(input == '-'){
                scale = std::max(scale/2,1);
                redraw = true;
                break;
            }
            if(input == 'w'){
                scrolHight = std::min(scrolHight+25,1000);
            }
            if(input == 'd'){
                scrolWidth = std::min(scrolWidth+25,1000);
            }
            if(input == 's'){
                scrolHight = std::max(scrolHight-25,0);
            }
            if(input == 'a'){
                scrolWidth = std::max(scrolWidth-25,0);
            }
            if(input == 'j'){
                if(labeling[labelingIdx] != _labelstart) {
                    labeling[labelingIdx] = (Label) (labeling[labelingIdx] - 1);

                    for(auto l : labeling){
                        std::cout << l << " ";
                    }
                    std::cout << std::endl;

                    redraw = true;
                    break;
                }
            }
            if(input == 'l'){
                if(labeling[labelingIdx] != _labelend) {
                    labeling[labelingIdx] = (Label) (labeling[labelingIdx] + 1);

                    for(auto l : labeling){
                        std::cout << l << " ";
                    }
                    std::cout << std::endl;

                    redraw = true;
                    break;
                }
            }
            if(input == 'i'){
                if(labelingIdx > 0){
                    labelingIdx -= 1;
                }
            }
            if(input == 'k'){
                if(labelingIdx < labeling.size()-1){
                    labelingIdx += 1;
                }
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}
