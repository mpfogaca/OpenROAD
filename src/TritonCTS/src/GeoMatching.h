#ifndef GEO_MATCHING
#define GEO_MATCHING

#include "Util.h"

#include <vector>
#include <limits>

namespace TritonCTS {

class Matching {
public:
        Matching(unsigned p0, unsigned p1) : _p0(p0), _p1(p1) {}

        unsigned getP0() const { return _p0; }
        unsigned getP1() const { return _p1; }
private:
        unsigned _p0;
        unsigned _p1;
};

class GeoMatching {
public:
        GeoMatching() = default;

        void addPoint(double x, double y) {  _points.emplace_back(x, y); }
        void run();
        unsigned getNumPoints() const { return _points.size(); }
        
        const std::vector<Matching>& allMatchings() const { return _matchings; } 

private:
        void normalizePoints();
        void computeAllThetas();
        void sortPoints();
        void findBestMatching();
        void writePlotFile();
        
        double computeTheta(double x, double y) const;
        unsigned numVertex(unsigned x, unsigned y) const;
        
        bool isOne(double pos) const { return (1 - pos) < std::numeric_limits<double>::epsilon(); }
        bool isZero(double pos) const { return pos < std::numeric_limits<double>::epsilon(); }
        
        std::vector<Point<double>> _points;
        std::vector<std::pair<double, unsigned>> _thetaIndexVector;
        std::vector<Matching> _matchings;
};

}

#endif 
