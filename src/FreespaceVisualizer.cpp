
#include "FreespaceVisualizer.h"

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

using std::cout;
using std::endl;
using namespace cv;
/*
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

        for (int y = 0; y < freespace.ySize(); ++y) {
            for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                auto cellwrapper = freespace.cell(y, xidx);
                int x = cellwrapper->x;
                //auto cell = cellwrapper->data;
                if (cellwrapper->data->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //connected cells lines

                Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                Point2d br = origin + Point2d(1.0, 0.0) * CS;

                rectangle(img, br, tl, Scalar(172, 172, 172), -1);
            }
        }

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

        for (int y = 0; y < freespace.ySize(); ++y) {
            for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                auto cellwrapper = freespace.cell(y, xidx);
                int x = cellwrapper->x;
                //auto cell = cellwrapper->data;
                if (cellwrapper->data->is_empty())
                    continue;

                //stupid mapping only flip y
                int imageY = (freespace.ySize() + 1 - y) * CS;

                Point2d origin = Point2d((x + 1) * CS, imageY);

                //connected cells lines

                Point2d tl = origin + Point2d(0.0,-1.0)*CS;
                Point2d tr = origin + Point2d(1.0,-1.0)*CS;
                Point2d bl = origin + Point2d(0.0,0.0)*CS;
                Point2d br = origin + Point2d(1.0,0.0)*CS;

                if(cellwrapper->upId != -1){
                    line(img, tl, tr, Scalar(255, 0, 0), 4);
                }
                if(cellwrapper->downId != -1){
                    line(img, bl, br, Scalar(255, 0, 0), 4);
                }
                if(cellwrapper->rightId != -1){
                    line(img, tr, br, Scalar(255, 0, 0), 4);
                }
                if(cellwrapper->leftId != -1){
                    line(img, tl, bl, Scalar(255, 0, 0), 4);
                }

                //eight points, eight lines

                Point2d p1 = Point2d(cellwrapper->data->leftPair.first.x, -cellwrapper->data->leftPair.first.y) * CS + origin;
                Point2d p2 = Point2d(cellwrapper->data->leftPair.second.x, -cellwrapper->data->leftPair.second.y) * CS + origin;
                Point2d p3 = Point2d(cellwrapper->data->topPair.first.x, -cellwrapper->data->topPair.first.y) * CS + origin;
                Point2d p4 = Point2d(cellwrapper->data->topPair.second.x, -cellwrapper->data->topPair.second.y) * CS + origin;
                Point2d p6 = Point2d(cellwrapper->data->rightPair.first.x, -cellwrapper->data->rightPair.first.y) * CS + origin;
                Point2d p5 = Point2d(cellwrapper->data->rightPair.second.x, -cellwrapper->data->rightPair.second.y) * CS + origin;
                Point2d p8 = Point2d(cellwrapper->data->bottomPair.first.x, -cellwrapper->data->bottomPair.first.y) * CS + origin;
                Point2d p7 = Point2d(cellwrapper->data->bottomPair.second.x, -cellwrapper->data->bottomPair.second.y) * CS + origin;

                line(img, p1, p2, Scalar(0, 0, 255), 2);
                line(img, p2, p3, Scalar(0, 0, 255), 2);
                line(img, p3, p4, Scalar(0, 0, 255), 2);
                line(img, p4, p5, Scalar(0, 0, 255), 2);
                line(img, p5, p6, Scalar(0, 0, 255), 2);
                line(img, p6, p7, Scalar(0, 0, 255), 2);
                line(img, p7, p8, Scalar(0, 0, 255), 2);
                line(img, p8, p1, Scalar(0, 0, 255), 2);


                //for (auto p: cell->importantUpYs) {
                //    Point2d ip = Point2d(p.first.x, -p.first.y) * CS + origin;
                //    circle(img, ip, CR, Scalar(0, 255, 0), -1);
                //}
                //for (auto p: cell->importantDownYs) {
                //    Point2d ip = Point2d(p.first.x, -p.first.y) * CS + origin;
                //    circle(img, ip, CR, Scalar(255, 0, 0), -1);
                //}

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
*/
int getRealY(int y, int maxY){
    return (maxY - y);
}
/*
void FreeSpacesVisualizer::show() {
//this will be awful
    int CS = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<int> partialYSums;
    int ySize = 1;
    for(auto& col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
    }

    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto& fs : freespaces[0]){
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
        int num = 0;

        Mat img(xSize*CS, ySize*CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                num++;
                //std::cout << "drawing " << xi << " " << y1 << std::endl;
                //freespaceorigin
                int xOrigin = partialXSums[xi];
                int yOrigin = partialYSums[yi];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code
                SparseFreespace &freespace = freespaces[yi][xi];

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                        Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                        Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                        Point2d br = origin + Point2d(1.0, 0.0) * CS;

                        rectangle(img, br, tl, Scalar(172, 172, 172), -1);
                    }
                }

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1)+fsOrigin, Point2d(x, -y2)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS)+fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0,-1.0)*CS;
                        Point2d tr = origin + Point2d(1.0,-1.0)*CS;
                        Point2d bl = origin + Point2d(0.0,0.0)*CS;
                        Point2d br = origin + Point2d(1.0,0.0)*CS;

                        if(cellwrapper->upId != -1){
                            line(img, tl, tr, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->downId != -1){
                            line(img, bl, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->rightId != -1){
                            line(img, tr, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->leftId != -1){
                            line(img, tl, bl, Scalar(255, 0, 0), 4);
                        }

                        //eight points, eight lines

                        Point2d p1 = Point2d(cellwrapper->data->leftPair.first.x, -cellwrapper->data->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cellwrapper->data->leftPair.second.x, -cellwrapper->data->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cellwrapper->data->topPair.first.x, -cellwrapper->data->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cellwrapper->data->topPair.second.x, -cellwrapper->data->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cellwrapper->data->rightPair.first.x, -cellwrapper->data->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cellwrapper->data->rightPair.second.x, -cellwrapper->data->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cellwrapper->data->bottomPair.first.x, -cellwrapper->data->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cellwrapper->data->bottomPair.second.x, -cellwrapper->data->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 2);
                        line(img, p2, p3, Scalar(0, 0, 255), 2);
                        line(img, p3, p4, Scalar(0, 0, 255), 2);
                        line(img, p4, p5, Scalar(0, 0, 255), 2);
                        line(img, p5, p6, Scalar(0, 0, 255), 2);
                        line(img, p6, p7, Scalar(0, 0, 255), 2);
                        line(img, p7, p8, Scalar(0, 0, 255), 2);
                        line(img, p8, p1, Scalar(0, 0, 255), 2);


                    }
                }
            }
        }
        std::cout << "drew all " << num << " populated freespaces of " << freespaces.size()*freespaces.size() << " possible\n";

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
 */
/*
void FreeSpacesVisualizer::showCandidates(std::vector<Candidate> candidates) {
    //this will be awful
    int CS = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<int> partialYSums;
    int ySize = 1;
    for(auto& col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
    }

    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto& fs : freespaces[0]){
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
                std::cout << "drawing " << xi << " " << y1 << std::endl;
                //freespaceorigin
                int xOrigin = partialXSums[xi];
                int yOrigin = partialYSums[yi];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code
                SparseFreespace &freespace = freespaces[yi][xi];

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                        Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                        Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                        Point2d br = origin + Point2d(1.0, 0.0) * CS;

                        rectangle(img, br, tl, Scalar(172, 172, 172), -1);
                    }
                }

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1)+fsOrigin, Point2d(x, -y2)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS)+fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0,-1.0)*CS;
                        Point2d tr = origin + Point2d(1.0,-1.0)*CS;
                        Point2d bl = origin + Point2d(0.0,0.0)*CS;
                        Point2d br = origin + Point2d(1.0,0.0)*CS;

                        if(cellwrapper->upId != -1){
                            line(img, tl, tr, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->downId != -1){
                            line(img, bl, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->rightId != -1){
                            line(img, tr, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->leftId != -1){
                            line(img, tl, bl, Scalar(255, 0, 0), 4);
                        }

                        //eight points, eight lines

                        Point2d p1 = Point2d(cellwrapper->data->leftPair.first.x, -cellwrapper->data->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cellwrapper->data->leftPair.second.x, -cellwrapper->data->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cellwrapper->data->topPair.first.x, -cellwrapper->data->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cellwrapper->data->topPair.second.x, -cellwrapper->data->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cellwrapper->data->rightPair.first.x, -cellwrapper->data->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cellwrapper->data->rightPair.second.x, -cellwrapper->data->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cellwrapper->data->bottomPair.first.x, -cellwrapper->data->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cellwrapper->data->bottomPair.second.x, -cellwrapper->data->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 2);
                        line(img, p2, p3, Scalar(0, 0, 255), 2);
                        line(img, p3, p4, Scalar(0, 0, 255), 2);
                        line(img, p4, p5, Scalar(0, 0, 255), 2);
                        line(img, p5, p6, Scalar(0, 0, 255), 2);
                        line(img, p6, p7, Scalar(0, 0, 255), 2);
                        line(img, p7, p8, Scalar(0, 0, 255), 2);
                        line(img, p8, p1, Scalar(0, 0, 255), 2);

                    }
                }
            }
        }

        //candidates

        for(int i=0;i<candidates.size();++i){
            auto candidate = candidates[i];
            bool last = i == candidates.size()-1;
            double y1 = ySize*CS - ((double)(candidate.getBegin().getPoint()) + candidate.getBegin().getFraction() + partialYSums[candidate.getCurveIndex()]+1)*CS;
            double y2 = ySize*CS - ((double)(candidate.getEnd().getPoint()) + candidate.getEnd().getFraction() + partialYSums[candidate.getCurveIndex()]+1)*CS;
            line(img, Point2d(0,y1), Point2d(xSize*CS,y1), Scalar(255, 0, 0), 1);
            line(img, Point2d(0,y2), Point2d(xSize*CS,y2), Scalar(255, 0, 0), 1);

            for(auto cov : candidate.visualMatching){
                int targetCurveIdx = cov.getCurveIndex();
                double x1 = (cov.getBegin().getFraction() + cov.getBegin().getPoint() + partialXSums[targetCurveIdx]+1)*CS;
                double x2 = (cov.getEnd().getFraction() + cov.getEnd().getPoint() + partialXSums[targetCurveIdx]+1)*CS;

                line(img, Point2d(x1,y1), Point2d(x2,y2), Scalar(0, last?255:0, last?0:255), 3);

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

        case 1:
            return {0,0,0};
        case 2:
            return {0,255,255};
        case 3:
            return {64,190,255};
        case 4:
            return {128,128,255};
        case 5:
            return {190,64,255};
        case 6:
            return {255,0,255};
        case 7:
            return {255,64,190};
        case 8:
            return {255,128,128};
        case 9:
            return {255,190,64};
        case 10:
            return {255,255,0};
        case 11:
            return {190,255,64};
        case 12:
            return {128,255,128};
        case 13:
            return {64,255,190};
        default:
            return {0,0,0};
    }
}

//TODO: not so fucking hacky
void
ClusteringVisulaizer::showClusteringStretched(Curves c, std::vector<std::vector<std::pair<Label, CPoint>>> groundthruth,
                                              std::vector<Candidate> candidates) {
    int scale = 8;

    int scrolHight = 0;
    int scrolWidth = 0;

    std::vector<Label> labeling(candidates.size());
    if(withAutocoloring){
        for(int i=0;i<candidates.size();++i){
            //tally up
            std::vector<int> labelcount(14,0);
            Candidate candidate = candidates[i];
            for(auto m : candidate.matching){
                for(CPoint s = m.getBegin();s<= m.getEnd();s.setPoint(s.getPoint()+1)){
                    int currentLabelIdx = 0;
                    while(groundthruth[m.getCurveIndex()][currentLabelIdx].second < s){
                        currentLabelIdx += 1;
                    }
                    labelcount[groundthruth[m.getCurveIndex()][currentLabelIdx].first] += 1;
                }
            }
            //find maximum
            Label maxLabel = 0;
            int count = -1;
            int totalcount = 0;
            for(int j=0;j<labelcount.size();j++){
                totalcount += labelcount[j];
                if(labelcount[j] > count){
                    count = labelcount[j];
                    maxLabel = (Label)j;
                }
            }
            if(count*2 > totalcount){
                labeling[i] = maxLabel;
            }
        }
    }
    int labelingIdx = 0;

    std::vector<int> partialLengths;
    partialLengths.push_back(1);
    for (auto curve: c) {
        partialLengths.push_back(partialLengths[partialLengths.size() - 1] + 1 + curve.size());
    }
    int xSize = (partialLengths[partialLengths.size() - 1]);
    int ySize = 1 + (14 - 0 - 1) + 5 + candidates.size() + 20 + 20 + 10;

    std::vector<std::vector<Label>> combinator;
    for(auto curve : c){
        combinator.emplace_back(curve.size()-1,0);
    }

    std::cout << xSize << " x " << ySize << std::endl;

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);

    double targetXSize = 200;
    double xScale = targetXSize/xSize;

    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(scale * ySize, scale * xSize*xScale, CV_8UC3, Scalar(255, 255, 255));

        for (int i = 0; i < c.size(); i++) {
            line(img, Point2d(partialLengths[i]*xScale, 1) * scale, Point2d((partialLengths[i + 1]-1)*xScale, 1) * scale, Scalar(0, 0, 0),
                 5);
            circle(img, Point2d(partialLengths[i]*xScale, 1) * scale, 5, Scalar(0, 0, 0), -1);
            circle(img, Point2d((partialLengths[i + 1]-1)*xScale, 1) * scale, 5, Scalar(0, 0, 0), -1);
        }

        for(Label label = (Label)(0+1);label!=14;label = (Label)(label+1)){
            double y = 1+label;
            line(img, Point2d(partialLengths[0]*xScale, y) * scale, Point2d(partialLengths[partialLengths.size()-1]*xScale, y) * scale, Scalar(0, 0, 0),
                 1);
        }

        for (int curveindex = 0; curveindex < c.size(); curveindex++) {
            CPoint start = {0, 0};
            for (auto assignment: groundthruth[curveindex]) {
                double y = 1 + assignment.first;

                double intxs = partialLengths[curveindex] + start.getPoint() + start.getFraction();
                double intxt = partialLengths[curveindex] + assignment.second.getPoint() + assignment.second.getFraction();

                rectangle(img,Point2d(intxs*xScale, ySize-32)*scale,Point2d(intxt*xScale, ySize-48)*scale, _labelColor(assignment.first),-1);
                rectangle(img,Point2d(intxs*xScale, ySize-32)*scale,Point2d(intxt*xScale, ySize-48)*scale, _labelColor(1),2);
                line(img, Point2d(intxs*xScale, y) * scale, Point2d(intxt*xScale, y) * scale, _labelColor(assignment.first), 5);
                circle(img, Point2d(intxs*xScale, y) * scale, 5, _labelColor(assignment.first), -1);
                circle(img, Point2d(intxt*xScale, y) * scale, 5, _labelColor(assignment.first), -1);

                start = assignment.second;
            }
        }

        //clear combinator
        combinator.clear();
        for(auto curve : c){
            combinator.emplace_back(curve.size()-1,0);
        }


        for (int i = 0; i < candidates.size(); i++) {
            auto wrappercandidate = candidates[i];
            Candidate candidate = wrappercandidate;
            double y = 1 + (14 - 0 - 1) + 5 + i;
            line(img, Point2d(partialLengths[0]*xScale, y) * scale, Point2d(partialLengths[partialLengths.size()-1]*xScale, y) * scale, Scalar(0, 0, 0),
                 1);
            for (auto matching: candidate.visualMatching) {

                for(int idx = matching.getBegin().getPoint();idx<=matching.getEnd().getPoint();idx++){
                    if(labeling[i] == 0){
                        continue;
                    }
                    if(combinator[matching.getCurveIndex()][idx] == 0 || combinator[matching.getCurveIndex()][idx] == labeling[i]){
                        combinator[matching.getCurveIndex()][idx] = labeling[i];
                    }else{
                        combinator[matching.getCurveIndex()][idx] = 1;
                    }
                }

                double intxs = partialLengths[matching.getCurveIndex()] + matching.getBegin().getFraction() + matching.getBegin().getPoint();
                double intxt = partialLengths[matching.getCurveIndex()] + matching.getEnd().getFraction() + matching.getEnd().getPoint();

                line(img, Point2d(intxs*xScale, y+0.1) * scale, Point2d(intxt*xScale, y-0.1) * scale, _labelColor(labeling[i]), 5);
                circle(img, Point2d(intxs*xScale, y+0.1) * scale, 5, _labelColor(labeling[i]), -1);
                circle(img, Point2d(intxt*xScale, y-0.1) * scale, 5, _labelColor(labeling[i]), -1);
            }
        }

        //draw combinator
        for (int i = 0; i < c.size(); i++) {
            int startx = partialLengths[i];
            for (int j=0;j<combinator[i].size();j++){
                //draw labeling of vertex
                rectangle(img,Point2d((startx+j)*xScale,ySize-12)*scale,Point2d((startx+j+1)*xScale,ySize-28)*scale, _labelColor(combinator[i][j]),-1);
            }
            rectangle(img,Point2d(startx*xScale,ySize-12)*scale,Point2d((startx+combinator[i].size())*xScale,ySize-28)*scale, _labelColor(1),2);
            //line(img, Point2d(partialLengths[i], 1) * scale, Point2d(partialLengths[i + 1]-1, 1) * scale, Scalar(0, 0, 0),5);
        }

        for (int i = 0; i < candidates.size(); i++) {
            auto wrappercandidate = candidates[i];
            Candidate candidate = wrappercandidate;
            for (auto matching: candidate.matching) {

                double intxs = partialLengths[matching.getCurveIndex()] + matching.getBegin().getPoint();// + matching.start.t;
                double intxt = partialLengths[matching.getCurveIndex()] + matching.getEnd().getPoint();// + matching.end.t;

                if(labeling[i]!=0) {
                    line(img, Point2d(intxs * xScale, ySize - 28) * scale, Point2d(intxs * xScale, ySize - 12) * scale,
                         _labelColor(1), 2);
                    //line(img, Point2d(intxt*xScale, ySize - 28) * scale, Point2d(intxt*xScale, ySize - 12) * scale,
                    //     _labelColor(transition), 2);
                }
                //circle(img, Point2d(intxs, y+0.1) * scale, 5, _labelColor(labeling[i]), -1);
                //circle(img, Point2d(intxt, y-0.1) * scale, 5, _labelColor(labeling[i]), -1);
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
                if(labeling[labelingIdx] != 0) {
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
                if(labeling[labelingIdx] != 14) {
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
*/
void SparseFreeSpacesVisualizer::show() {
//this will be awful
    double CS = 32.0;
    while(CS/8.0 > 1.0/(double)(freespaces.size())){
        CS /= 2.0;
    }

    std::cout << std::flush;
    int scrolHight = 0;
    int scrolWidth = 0;

    //assumes its a square
    std::vector<int> partialYSums;
    int ySize = 1;
    std::vector<int> partialXSums;
    int xSize = 1;
    for(auto& col : freespaces){
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
        partialXSums.push_back(xSize);
        xSize += (col.front().ySize() + 2);
    }

    namedWindow("+/-: zoom, w/a/s/d: move, q: quit", WINDOW_NORMAL);
    //namedWindow("controlWin", WINDOW_AUTOSIZE);

    //createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    //createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while(redraw) {
        redraw = false;

        Mat img(xSize*CS, ySize*CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                SparseFreespace &freespace = freespaces[yi][xi];
                int actualX = freespace.TID;
                int actualY = freespace.BID;
                //std::cout << "drawing " << actualX << " " << actualY << std::endl;
                //freespaceorigin
                int xOrigin = partialXSums[actualX];
                int yOrigin = partialYSums[actualY];


                Point2d fsOrigin = Point2d(xOrigin*CS, getRealY(yOrigin, ySize)*CS);
                //circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS) + fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                        Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                        Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                        Point2d br = origin + Point2d(1.0, 0.0) * CS;

                        rectangle(img, br, tl, Scalar(172, 172, 172), -1);
                    }
                }

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1)+fsOrigin, Point2d(x, -y2)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y)+fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y-1)*CS)+fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0,-1.0)*CS;
                        Point2d tr = origin + Point2d(1.0,-1.0)*CS;
                        Point2d bl = origin + Point2d(0.0,0.0)*CS;
                        Point2d br = origin + Point2d(1.0,0.0)*CS;

                        if(cellwrapper->upId != -1){
                            line(img, tl, tr, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->downId != -1){
                            line(img, bl, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->rightId != -1){
                            line(img, tr, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->leftId != -1){
                            line(img, tl, bl, Scalar(255, 0, 0), 4);
                        }

                        //eight points, eight lines

                        Point2d p1 = Point2d(cellwrapper->data->leftPair.first.x, -cellwrapper->data->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cellwrapper->data->leftPair.second.x, -cellwrapper->data->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cellwrapper->data->topPair.first.x, -cellwrapper->data->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cellwrapper->data->topPair.second.x, -cellwrapper->data->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cellwrapper->data->rightPair.first.x, -cellwrapper->data->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cellwrapper->data->rightPair.second.x, -cellwrapper->data->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cellwrapper->data->bottomPair.first.x, -cellwrapper->data->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cellwrapper->data->bottomPair.second.x, -cellwrapper->data->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 2);
                        line(img, p2, p3, Scalar(0, 0, 255), 2);
                        line(img, p3, p4, Scalar(0, 0, 255), 2);
                        line(img, p4, p5, Scalar(0, 0, 255), 2);
                        line(img, p5, p6, Scalar(0, 0, 255), 2);
                        line(img, p6, p7, Scalar(0, 0, 255), 2);
                        line(img, p7, p8, Scalar(0, 0, 255), 2);
                        line(img, p8, p1, Scalar(0, 0, 255), 2);
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
            imshow("+/-: zoom, w/a/s/d: move, q: quit", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if(input == '+'){
                CS = std::min(CS*2,128.0);
                redraw = true;
                break;
            }
            if(input == '-'){
                CS = std::max(CS/2,1.0/1024);
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
/*
void SparseFreeSpacesVisualizer::showCandidates(std::vector<Candidate> candidates) {
    std::cout << "test";
    return;
//this will be awful
    double CS = 32.0;
    while(CS/8.0 > 1.0/(double)(freespaces.size())){
        CS /= 2.0;
    }

    int scrolHight = 0;
    int scrolWidth = 0;

    //assumes its a square
    std::vector<int> partialYSums;
    int ySize = 1;
    std::vector<int> partialXSums;
    int xSize = 1;
    for (auto& col: freespaces) {
        partialYSums.push_back(ySize);
        ySize += (col.front().ySize() + 2);
        partialXSums.push_back(xSize);
        xSize += (col.front().ySize() + 2);
    }

    namedWindow("winImage", WINDOW_NORMAL);
    namedWindow("controlWin", WINDOW_AUTOSIZE);

    createTrackbar("Hscroll", "controlWin", &scrolHight, 1000);
    createTrackbar("Wscroll", "controlWin", &scrolWidth, 1000);


    bool redraw = true;
    while (redraw) {
        redraw = false;

        Mat img(xSize * CS, ySize * CS, CV_8UC3, Scalar(255, 255, 255));

        for (int yi = 0; yi < freespaces.size(); yi++) {
            for (int xi = 0; xi < freespaces[yi].size(); xi++) {
                SparseFreespace &freespace = freespaces[yi][xi];
                int actualX = freespace.TID;
                int actualY = freespace.BID;
                std::cout << "drawing " << actualX << " " << actualY << std::endl;
                //freespaceorigin
                int xOrigin = partialXSums[actualX];
                int yOrigin = partialYSums[actualY];


                Point2d fsOrigin = Point2d(xOrigin * CS, getRealY(yOrigin, ySize) * CS);
                //circle(img,fsOrigin, 5, Scalar(0, 0, 0), -1);


                //now copy code

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y - 1) * CS) + fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                        Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                        Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                        Point2d br = origin + Point2d(1.0, 0.0) * CS;

                        rectangle(img, br, tl, Scalar(172, 172, 172), -1);
                    }
                }

                //draw vertexlines
                for (int i = 0; i <= freespace.xSize(); ++i) {
                    int x = (i + 1) * CS;
                    int y1 = CS;
                    int y2 = (freespace.ySize() + 1) * CS;
                    line(img, Point2d(x, -y1) + fsOrigin, Point2d(x, -y2) + fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int i = 0; i <= freespace.ySize(); ++i) {
                    int y = (i + 1) * CS;
                    int x1 = CS;
                    int x2 = (freespace.xSize() + 1) * CS;
                    line(img, Point2d(x1, -y) + fsOrigin, Point2d(x2, -y) + fsOrigin, Scalar(0, 0, 0), 3);
                }

                for (int y = 0; y < freespace.ySize(); ++y) {
                    for (int xidx = 0; xidx < freespace.row(y).size(); ++xidx) {

                        auto cellwrapper = freespace.cell(y, xidx);
                        int x = cellwrapper->x;
                        //auto cell = cellwrapper->data;
                        if (cellwrapper->data->is_empty())
                            continue;

                        //stupid mapping only flip y
                        //int imageY = (freespace.ySize() + 1 - y) * CS;

                        Point2d origin = Point2d((x + 1) * CS, (-y - 1) * CS) + fsOrigin;

                        //connected cells lines

                        Point2d tl = origin + Point2d(0.0, -1.0) * CS;
                        Point2d tr = origin + Point2d(1.0, -1.0) * CS;
                        Point2d bl = origin + Point2d(0.0, 0.0) * CS;
                        Point2d br = origin + Point2d(1.0, 0.0) * CS;

                        if(cellwrapper->upId != -1){
                            line(img, tl, tr, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->downId != -1){
                            line(img, bl, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->rightId != -1){
                            line(img, tr, br, Scalar(255, 0, 0), 4);
                        }
                        if(cellwrapper->leftId != -1){
                            line(img, tl, bl, Scalar(255, 0, 0), 4);
                        }

                        //eight points, eight lines

                        Point2d p1 = Point2d(cellwrapper->data->leftPair.first.x, -cellwrapper->data->leftPair.first.y) * CS + origin;
                        Point2d p2 = Point2d(cellwrapper->data->leftPair.second.x, -cellwrapper->data->leftPair.second.y) * CS + origin;
                        Point2d p3 = Point2d(cellwrapper->data->topPair.first.x, -cellwrapper->data->topPair.first.y) * CS + origin;
                        Point2d p4 = Point2d(cellwrapper->data->topPair.second.x, -cellwrapper->data->topPair.second.y) * CS + origin;
                        Point2d p6 = Point2d(cellwrapper->data->rightPair.first.x, -cellwrapper->data->rightPair.first.y) * CS + origin;
                        Point2d p5 = Point2d(cellwrapper->data->rightPair.second.x, -cellwrapper->data->rightPair.second.y) * CS + origin;
                        Point2d p8 = Point2d(cellwrapper->data->bottomPair.first.x, -cellwrapper->data->bottomPair.first.y) * CS + origin;
                        Point2d p7 = Point2d(cellwrapper->data->bottomPair.second.x, -cellwrapper->data->bottomPair.second.y) * CS + origin;

                        line(img, p1, p2, Scalar(0, 0, 255), 2);
                        line(img, p2, p3, Scalar(0, 0, 255), 2);
                        line(img, p3, p4, Scalar(0, 0, 255), 2);
                        line(img, p4, p5, Scalar(0, 0, 255), 2);
                        line(img, p5, p6, Scalar(0, 0, 255), 2);
                        line(img, p6, p7, Scalar(0, 0, 255), 2);
                        line(img, p7, p8, Scalar(0, 0, 255), 2);
                        line(img, p8, p1, Scalar(0, 0, 255), 2);


                    }
                }
            }
        }

        //candidates

        for(int i=0;i<candidates.size();++i){
            auto candidate = candidates[i];
            bool last = i == candidates.size()-1;
            double y1 = ySize*CS - ((double)(candidate.getBegin().getPoint()) + candidate.getBegin().getFraction() + partialYSums[candidate.getCurveIndex()]+1)*CS;
            double y2 = ySize*CS - ((double)(candidate.getEnd().getPoint()) + candidate.getEnd().getFraction() + partialYSums[candidate.getCurveIndex()]+1)*CS;
            line(img, Point2d(0,y1), Point2d(xSize*CS,y1), Scalar(255, 0, 0), 1);
            line(img, Point2d(0,y2), Point2d(xSize*CS,y2), Scalar(255, 0, 0), 1);

            for(auto cov : candidate.visualMatching){
                int targetCurveIdx = cov.getCurveIndex();
                double x1 = (cov.getBegin().getFraction() + cov.getBegin().getPoint() + partialXSums[targetCurveIdx]+1)*CS;
                double x2 = (cov.getEnd().getFraction() + cov.getEnd().getPoint() + partialXSums[targetCurveIdx]+1)*CS;

                line(img, Point2d(x1,y1), Point2d(x2,y2), Scalar(0, last?255:0, last?0:255), 3);

            }


        }

        int winH = 1800;
        int winW = 3200;
        if (winH >= img.rows)winH = img.rows - 1;
        if (winW >= img.cols)winW = img.cols - 1;
        while (true) {
            int truescrolHight = (img.rows - winH) * scrolHight / 1000;
            int truescrolWidth = (img.cols - winW) * scrolWidth / 1000;
            Mat winImage = img(Rect(truescrolWidth, img.rows - winH - truescrolHight, winW, winH));
            imshow("winImage", winImage);
            int input = waitKey(0);
            if (input == 'q')
                break;
            if (input == '+') {
                CS = std::min(CS * 2, 128.0);
                redraw = true;
                break;
            }
            if (input == '-') {
                CS = std::max(CS / 2, 1.0 / 1024);
                redraw = true;
                break;
            }
            if (input == 'w') {
                scrolHight = std::min(scrolHight + 25, 1000);
            }
            if (input == 'd') {
                scrolWidth = std::min(scrolWidth + 25, 1000);
            }
            if (input == 's') {
                scrolHight = std::max(scrolHight - 25, 0);
            }
            if (input == 'a') {
                scrolWidth = std::max(scrolWidth - 25, 0);
            }
        }
    }

    cv::destroyAllWindows();

    cv::waitKey(1);
}
*/