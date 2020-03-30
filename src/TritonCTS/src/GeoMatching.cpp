#include "GeoMatching.h"
#include <tuple>
#include <iostream>
#include <cmath> 
#include <algorithm>
#include <fstream>

namespace TritonCTS {

void GeoMatching::normalizePoints() {
        double xMax = -std::numeric_limits<double>::infinity();
        double xMin =  std::numeric_limits<double>::infinity();
        double yMax = -std::numeric_limits<double>::infinity();
        double yMin =  std::numeric_limits<double>::infinity();
        for (const Point<double>& p : _points) {         
                xMax = std::max(p.getX(), xMax);        
                yMax = std::max(p.getY(), yMax);        
                xMin = std::min(p.getX(), xMin);        
                yMin = std::min(p.getY(), yMin);        
        }

        double xSpan = xMax - xMin;
        double ySpan = yMax - yMin;
        for (Point<double>& p : _points) {
                double x = p.getX();
                double xNorm = (x - xMin) / xSpan; 
                double y = p.getY();
                double yNorm = (y - yMin)/ ySpan; 
                p = Point<double>(xNorm, yNorm);
        }
}

void GeoMatching::computeAllThetas() {
        for (unsigned idx = 0; idx < _points.size(); ++idx) {
                const Point<double>& p = _points[idx];
                double theta = computeTheta(p.getX(), p.getY());
                _thetaIndexVector.emplace_back(theta, idx);
        }
}

void GeoMatching::sortPoints() {
        std::sort(_thetaIndexVector.begin(), _thetaIndexVector.end());
}

double GeoMatching::computeTheta(double x, double y) const {
        if (isOne(x) && isOne(y)) {
                return 0.5;        
        }
        
        unsigned quad = numVertex(std::min(unsigned(2.0 * x), (unsigned) 1), 
                                  std::min(unsigned(2.0 * y), (unsigned) 1));

        double t = computeTheta(2 * std::fabs(x - 0.5),
                                2 * std::fabs(y - 0.5)); 

        if (quad % 2 == 1) {
                t = 1 - t;
        }

        double integral;
        double fractal = std::modf((quad + t) / 4.0 + 7.0 / 8.0, &integral);
        return fractal;
}

unsigned GeoMatching::numVertex(unsigned x, unsigned y) const {
        if ((x == 0) && (y == 0)) {
                return 0;
        } else if ((x == 0) && (y == 1)) {
                return 1;
        } else if ((x == 1) && (y == 1)) {
                return 2; 
        } else if ((x == 1) && (y == 0)) {
                return 3;
        }
        
        std::cout << "[ERROR] Invalid parameters in " << __func__ << "\n";
        std::exit(1);
}

void GeoMatching::findBestMatching() {
        std::vector<Matching> matching0;
        std::vector<Matching> matching1;
        double matchingCost0 = 0.0;
        double matchingCost1 = 0.0;

        for (unsigned i = 0; i < _thetaIndexVector.size() - 1; ++i) {
                unsigned idx0 = _thetaIndexVector[i].second;
                unsigned idx1 = _thetaIndexVector[i + 1].second;
                Point<double>& p0 = _points[idx0];
                Point<double>& p1 = _points[idx1];
                
                double cost = p0.computeDist(p1);
                if (i % 2 == 0) {
                        matching0.emplace_back(idx0, idx1);
                        matchingCost0 += cost;
                } else {
                        matching1.emplace_back(idx0, idx1);
                        matchingCost1 += cost;
                }
        }
        
        // Cost 1 needs to sum the distance from first to last
        // element
        unsigned idx0 = _thetaIndexVector.front().second;
        unsigned idx1 = _thetaIndexVector.back().second; 
        matching1.emplace_back(idx0, idx1);
        Point<double>& p0 = _points[idx0];
        Point<double>& p1 = _points[idx1];
        matchingCost1 += p0.computeDist(p1);
        
        std::cout << "Matching0 size: " << matching0.size() << "\n";
        std::cout << "Matching1 size: " << matching1.size() << "\n";
        if (matchingCost0 < matchingCost1) {
                _matchings = matching0;
        } else {
                _matchings = matching1;
        }
}

void GeoMatching::run() {
        normalizePoints();
        computeAllThetas();
        sortPoints();
        findBestMatching();        
        writePlotFile();
        //exit(1);
}

void GeoMatching::writePlotFile() {
        std::ofstream file("plot.py");
        file << "import numpy as np\n"; 
        file << "import matplotlib.pyplot as plt\n";
        file << "import matplotlib.path as mpath\n";
        file << "import matplotlib.lines as mlines\n";
        file << "import matplotlib.patches as mpatches\n";
        file << "from matplotlib.collections import PatchCollection\n\n";

        for (unsigned idx = 0; idx < _thetaIndexVector.size() - 1; idx += 2) {
                unsigned idx0 = _thetaIndexVector[idx].second;
                unsigned idx1 = _thetaIndexVector[idx+1].second;
                Point<double>& p0 = _points[idx0];
                Point<double>& p1 = _points[idx1];
                file << "plt.scatter(" << p0.getX() << ", " << p0.getY() << ")\n";
                file << "plt.scatter(" << p1.getX() << ", " << p1.getY() << ")\n";
                file << "plt.plot([" << p0.getX() << ", " << p1.getX() << "], ["
                     << p0.getY() << ", " << p1.getY() << "])\n";
        }

        file << "plt.show()\n"; 
        file.close();
}

}
