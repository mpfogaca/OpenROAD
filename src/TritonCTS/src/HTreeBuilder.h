////////////////////////////////////////////////////////////////////////////////////
// Authors: Mateus Fogaca
//          (Ph.D. advisor: Ricardo Reis)
//          Jiajia Li
//          Andrew Kahng
// Based on:
//          K. Han, A. B. Kahng and J. Li, "Optimal Generalized H-Tree Topology and 
//          Buffering for High-Performance and Low-Power Clock Distribution", 
//          IEEE Trans. on CAD (2018), doi:10.1109/TCAD.2018.2889756.
//
//
// BSD 3-Clause License
//
// Copyright (c) 2018, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////////////

#ifndef HTREEBUILDER_H
#define HTREEBUILDER_H

#include "TreeBuilder.h"
#include "CtsOptions.h"
#include "Util.h"

#include <limits>
#include <cmath>

namespace TritonCTS {

class SegmentBuilder {
public:
        SegmentBuilder(const std::string instPrefix,
                       const std::string netPrefix,
                       Point<double> root, 
                       Point<double> target,
                       const std::vector<unsigned>& techCharWires,
                       Clock& clock,
                       Clock::SubNet& drivingSubNet,
                       TechChar& techChar,
                       unsigned techCharDistUnit) :
                       _instPrefix(instPrefix),
                       _netPrefix(netPrefix),                       
                       _root(root), 
                       _target(target), 
                       _techCharWires(techCharWires), 
                       _clock(&clock),
                       _drivingSubNet(&drivingSubNet),
                       _techChar(&techChar),
                       _techCharDistUnit(techCharDistUnit) {}

        void build();
        void forceBufferInSegment(std::string master);
        Clock::SubNet* getDrivingSubNet() const { return _drivingSubNet; }

protected:
        const std::string     _instPrefix;       
        const std::string     _netPrefix;       
        Point<double>         _root;
        Point<double>         _target;
        std::vector<unsigned> _techCharWires;
        Clock*                _clock;
        Clock::SubNet*        _drivingSubNet;
        TechChar*             _techChar;
        unsigned              _techCharDistUnit;
        bool                  _forceBuffer;
        unsigned              _numBuffers = 0;

        void buildVerticalConnection();
        void buildHorizontalConnection();
        void buildLShapeConnection();
};

//-----------------------------------------------------------------------------

class TopologyWire {
public:
        TopologyWire(unsigned length) : _length(length) {};
        
        void addCharWire(unsigned key) { _charWires.push_back(key); }
        unsigned getLength() const { return _length; }
        void setLength(unsigned length) { _length = length; }
        unsigned getOutputSlew() const { return _outputSlew; }
        void setOutputSlew(unsigned slew) { _outputSlew = slew; }
        unsigned getOutputCap() const { return _outputCap; }
        void setOutputCap(unsigned cap) { _outputCap = cap; }
        const std::vector<unsigned>& allCharWires() const { return _charWires; }
protected:
        unsigned _length;
        unsigned _outputSlew;
        unsigned _outputCap;        
        std::vector<unsigned> _charWires;
};

//-----------------------------------------------------------------------------

class LevelTopology {
public:
        enum Direction : uint8_t {HORIZONTAL, VERTICAL};
        static constexpr unsigned NO_PARENT = std::numeric_limits<unsigned>::max();
        
        LevelTopology() = default;        
        LevelTopology(double length) : _length(length), _outputSlew(0), _outputCap(0) {};
               
        void addWireSegment(unsigned idx) { _wireSegments.push_back(idx); }
               
        unsigned addBranchingPoint(const Point<double>& loc, unsigned parent) { 
                _branchPointLoc.push_back(loc);
                _parents.push_back(parent);
                _branchSinkLocs.resize(_branchPointLoc.size());
                _branchDrivingSubNet.resize(_branchPointLoc.size(), nullptr);
                return _branchPointLoc.size() - 1;
        }

        void addSinkToBranch(unsigned branchIdx, const Point<double>& sinkLoc) {
                _branchSinkLocs[branchIdx].push_back(sinkLoc);                        
        }
 
        Point<double>& getBranchingPoint(unsigned idx) { 
                return _branchPointLoc[idx]; 
        }
 
        unsigned getBranchingPointParentIdx(unsigned idx) const {
                return _parents[idx];
        }
 
        double getLength() const { return _length; } 
 
        void forEachBranchingPoint(std::function<void(unsigned, Point<double>)> func) const {
                for (unsigned idx = 0; idx < _branchPointLoc.size(); ++idx) {
                        func(idx, _branchPointLoc[idx]);
                }
        }
 
        Clock::SubNet* getBranchDrivingSubNet(unsigned idx) const { 
                return _branchDrivingSubNet[idx]; 
        }
 
        void setBranchDrivingSubNet(unsigned idx, Clock::SubNet& subNet) {
                _branchDrivingSubNet[idx] = &subNet;
        }
 
        const std::vector<unsigned>& getWireSegments() const { return _wireSegments; }
        
        const std::vector<Point<double>>& getBranchSinksLocations(unsigned branchIdx) const {
               return _branchSinkLocs[branchIdx];
        }
 
        void setOutputSlew(unsigned slew) { _outputSlew = slew; }
        unsigned getOutputSlew() const { return _outputSlew; }
        void setOutputCap(unsigned cap) { _outputCap = cap; }
        unsigned getOutputCap() const { return _outputCap; }


        void setWidth(double width) { _width = width; }
        double getWidth() const { return _width; }
        void setHeight(double height) { _height = height; }
        double getHeight() const { return _height; }
        void setNumSinks(unsigned numSinks) { _numSinks = numSinks; }
        unsigned getNumSinks() const { return _numSinks; }
        void setBranchingFactor(unsigned factor) { _branchingFactor = factor; }
        unsigned getBranchingFactor() const { return _branchingFactor; }
        void setDirection(Direction dir) { _direction = dir; }
        bool isVertical() const { return _direction == VERTICAL; } 

        void createTopologyWire(unsigned length) { _topologyWires.push_back(length); }
        TopologyWire& getFirstWire() { return _topologyWires[0]; }
       
        void forEachTopologyWire(const std::function<void(unsigned, TopologyWire&)>& func) {
                for (unsigned idx = 0; idx < _topologyWires.size(); ++idx) {
                        func(idx, _topologyWires[idx]);
                } 
        }

        void forEachTopologyWire(const std::function<void(TopologyWire&)>& func) {
                for (TopologyWire& wire: _topologyWires) {
                        func(wire);
                } 
        }
        
        unsigned addBranchLocation(Point<double> loc) { 
                _branchLocs.push_back(loc); 
                unsigned numBranches = _branchLocs.size();
                unsigned idx = numBranches - 1;
                _branchSinkLocs.resize(numBranches);
                _parentBranchIdx.resize(numBranches);
                _topologyWireIdx.resize(numBranches);
                _branchDrivingSubNet.resize(numBranches, nullptr);
                return idx;
        }

        void setBranchLocation(unsigned idx, Point<double> loc) { _branchLocs[idx] = loc; } 
        Point<double> getBranchLocation(unsigned idx) const { return _branchLocs[idx]; }
        void setParentBranchIdx(unsigned idx, unsigned parent) { _parentBranchIdx[idx] = parent; }
        unsigned getParentBranchIdx(unsigned idx) const { return _parentBranchIdx[idx]; }
        void addUpstreamBranchIdx(unsigned idx) { _upstreamBranchIdx.push_back(idx); }
        unsigned getUpstreamBranchIdx(unsigned idx) const { return _upstreamBranchIdx[idx]; }
        void setTopologyWireIdx(unsigned branchIdx, unsigned wireIdx) { _topologyWireIdx[branchIdx] = wireIdx; }
        unsigned getTopologyWireIdx(unsigned branchIdx) const { return _topologyWireIdx[branchIdx]; }
        const TopologyWire& getTopologyWire(unsigned branchIdx) { 
                return _topologyWires[getTopologyWireIdx(branchIdx)]; 
        }
        void computeBranchLocations(LevelTopology* parentTopology);
        void refineBranchLocations(unsigned parentIdx,
                                   const std::vector<Point<double>>& sinkLocs,
                                   const Point<double> rootLoc);
        void forEachBranchLocation(std::function<void(unsigned, Point<double>)> func) const {
                for (unsigned idx = 0; idx < _branchLocs.size(); ++idx) {
                        func(idx, _branchLocs[idx]);
                }
        }
private:
        double                         _length;
        unsigned                       _outputSlew;
        unsigned                       _outputCap;
        std::vector<unsigned>          _wireSegments;
        std::vector<Point<double>>     _branchPointLoc;
        std::vector<unsigned>          _parents;

        double    _width;
        double    _height;
        unsigned  _numSinks;
        unsigned  _branchingFactor;
        Direction _direction;
        std::vector<Point<double>>  _branchLocs;
        std::vector<unsigned>       _upstreamBranchIdx;
        std::vector<unsigned>       _parentBranchIdx;
        std::vector<unsigned>       _topologyWireIdx;
        std::vector<TopologyWire>   _topologyWires;
        std::vector<Clock::SubNet*> _branchDrivingSubNet;
        std::vector<std::vector<Point<double>>> _branchSinkLocs;
};

//-----------------------------------------------------------------------------

class HTreeBuilder : public TreeBuilder {
public:
        HTreeBuilder(CtsOptions& options, Clock& net) : 
                     TreeBuilder(options, net) {};
        
        void run();

        void plotSolution();
private:
        void initSinkRegion();
        void initTopology();
        void reportTopology() const;
        unsigned computeNearestSegmentLength(double length);
        void computeCharWires();
        void computeCharWire(TopologyWire& wire, unsigned inputCap, unsigned inputSlew);
        void computeBranchLocations();

        void computeLevelTopology(unsigned level, double width, double height);
        unsigned computeNumberOfSinksPerSubRegion(unsigned level) const;
        void computeSubRegionSize(unsigned level, double& width, double& height) const;
        unsigned computeMinDelaySegment(unsigned length) const;
        unsigned computeMinDelaySegment(unsigned length, 
                                        unsigned inputSlew,
                                        unsigned inputCap,
                                        unsigned slewThreshold,
                                        unsigned tolerance,
                                        unsigned &outputSlew,
                                        unsigned &outputCap) const;
        void reportWireSegment(unsigned key) const;
        void createClockSubNets();
        void createSingleBufferClockNet();
        void initTopLevelSinks(std::vector<std::pair<float,float>>& sinkLocations);
        void computeBranchSinks(LevelTopology& topology, unsigned branchIdx,
                                std::vector<std::pair<float,float>>& sinkLocations) const;
        
        bool isVertical(unsigned level) const { return level % 2 == 0; }
        bool isHorizontal(unsigned level) const { return !isVertical(level); }

        unsigned computeGridSizeX(unsigned level) const { return std::pow(2, (level + 1) / 2); }
        unsigned computeGridSizeY(unsigned level) const { return std::pow(2, level / 2); }
        
        void computeBranchingPoints(unsigned level, LevelTopology& topology);
        void refineBranchingPointsWithClustering(LevelTopology& topology,
                                                 unsigned level,
                                                 unsigned branchPtIdx1,
                                                 unsigned branchPtIdx2, 
                                                 const Point<double>& rootLocation,
                                                 const std::vector<std::pair<float, float>>& topLevelSinks );
        
        bool isSubRegionTooSmall(double width, double height) const {
                if (width < _minLengthSinkRegion || height < _minLengthSinkRegion ) {
                        return true;
                }
                return false;       
        }

        bool isNumberOfSinksTooSmall(unsigned numSinksPerSubRegion) const {
                if (numSinksPerSubRegion < _numMaxLeafSinks) {
                        return true;
                }
                return false;
        }

protected:
        Box<double>                _sinkRegion;
        std::vector<LevelTopology> _topologyForEachLevel;
        
        DBU      _wireSegmentUnit     = -1;
        unsigned _minInputCap         =  1;
        unsigned _numMaxLeafSinks     =  0;
        unsigned _minLengthSinkRegion =  0;
        unsigned _clockTreeMaxDepth   =  0;
        unsigned _maxLevel            =  0;
}; 

}

#endif
