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

#include "HTreeBuilder.h"
#include "third_party/CKMeans/clustering.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

namespace TritonCTS {

void HTreeBuilder::initSinkRegion() {
        unsigned wireSegmentUnitInMicron = _techChar->getLengthUnit(); 
        DBU dbUnits = _options->getDbUnits();
        _wireSegmentUnit = wireSegmentUnitInMicron * dbUnits;

        std::cout << " Wire segment unit: " << _wireSegmentUnit << " dbu ("
                  << wireSegmentUnitInMicron << " um)\n";

        Box<DBU> sinkRegionDbu = _clock.computeSinkRegion();
        std::cout << " Original sink region: " << sinkRegionDbu << "\n";
        
        _sinkRegion = sinkRegionDbu.normalize(1.0/_wireSegmentUnit);
                std::cout << " Normalized sink region: " << _sinkRegion << "\n";
        std::cout << "    Width:  " << _sinkRegion.getWidth() << "\n";
        std::cout << "    Height: " << _sinkRegion.getHeight() << "\n";
}

void HTreeBuilder::run() {
        std::cout << " Generating H-Tree topology for net " << _clock.getName() << "...\n";
        std::cout << "    Tot. number of sinks: " << _clock.getNumSinks() << "\n";
       
        _minInputCap = _techChar->getActualMinInputCap();
        _numMaxLeafSinks = _options->getNumMaxLeafSinks();
        _minLengthSinkRegion = _techChar->getMinSegmentLength() * 2;

        initSinkRegion();
        _clockTreeDepth = computeClockTreeDepth();
        _topologyForEachLevel.resize(_clockTreeDepth);
        std::cout << " Clock tree depth: " << _clockTreeDepth << "\n";

        if(_clockTreeDepth < 1) {
                createSingleBufferClockNet(); 
                return;
        } 
        
        for (unsigned depth = _clockTreeDepth; depth >= 1; --depth) {
               computeLevelTopology(depth); 
        }

        std::cout << " Computing branch locations: Level ";
        for (unsigned depth = 1; depth <= _clockTreeDepth; ++depth) {
                computeBranchingPoints(depth);                        
        }
        std::cout << "\n";

        createClockSubNets();

}

unsigned HTreeBuilder::computeClockTreeDepth() {
        unsigned clockTreeMaxDepth = _options->getClockTreeMaxDepth();
        
        for (unsigned depth = 1; depth <= clockTreeMaxDepth; ++depth) {
                double width  = 0;
                double height = 0;
                computeSubRegionSize(depth, width, height);
                if (isSubRegionTooSmall(width, height)) {
                        if (_options->isFakeLutEntriesEnabled() && _minLengthSinkRegion) {
                                unsigned minIndex = 1;
                                _techChar->createFakeEntries(_minLengthSinkRegion, minIndex);
                                _minLengthSinkRegion = 2 * minIndex;
                        } else {
                                return depth - 1;
                        }
                }

                unsigned numSinksPerSubRegion = computeNumberOfSinksPerSubRegion(depth); 
                if (isNumberOfSinksTooSmall(numSinksPerSubRegion)) {
                        return depth;
                }         
        }          
        
        return 0;
}

inline 
unsigned HTreeBuilder::computeNumberOfSinksPerSubRegion(unsigned level) const {
        unsigned totalNumSinks = _clock.getNumSinks();
        unsigned numRoots = std::pow(2, level);
        double numSinksPerRoot = (double) totalNumSinks / numRoots;
        return (unsigned) std::ceil(numSinksPerRoot);
}

inline 
void HTreeBuilder::computeSubRegionSize(unsigned level, double& width, double& height) const {
        unsigned gridSizeX = computeGridSizeX(level);
        unsigned gridSizeY = computeGridSizeY(level);
        width = _sinkRegion.getWidth() / gridSizeX;
        height = _sinkRegion.getHeight() / gridSizeY;
}

void HTreeBuilder::computeLevelTopology(unsigned level) {
        unsigned numSinksPerSubRegion = computeNumberOfSinksPerSubRegion(level);
        
        double width  = 0.0;
        double height = 0.0;
        computeSubRegionSize(level, width, height);
        
        std::cout << " Level " << level << "\n";
        std::cout << "    Direction: " << ((isVertical(level)) ? ("Vertical") : ("Horizontal")) 
                  << "\n";
        std::cout << "    # sinks per sub-region: " << numSinksPerSubRegion << "\n";        
        std::cout << "    Sub-region size: " << width << " X " << height << "\n";      
      
        unsigned minLength = _minLengthSinkRegion / 2;
        unsigned segmentLength = std::round(width/(2.0*minLength))*minLength;
        if (isVertical(level)) {
                segmentLength = std::round(height/(2.0*minLength))*minLength; 
        }
        segmentLength = std::max<unsigned>(segmentLength, 1);        
        
        LevelTopology topology(segmentLength);
        std::cout << "    Segment length (rounded): " << segmentLength << "\n";
      
        unsigned load = 1;
        unsigned outSlew = 1; 
        if (level == _clockTreeDepth) {
                load = numSinksPerSubRegion;
                outSlew = _options->getMaxSlew();
        } else {
                load = _topologyForEachLevel[level].getInputCap();
                outSlew = _topologyForEachLevel[level].getInputSlew();                              
        }

        std::cout << "    Load: " << load; 
        std::cout << " Output slew: " << outSlew << "\n"; 
        
        unsigned currLength = 0;
        _techChar->forEachSegmentLengthReversed( [&] (unsigned charSegLength) {
                unsigned numWires = (segmentLength - currLength) / charSegLength;
                 
                if (numWires < 1) {
                        return;
                }
                
                currLength += numWires * charSegLength;
                for (unsigned wireCount = 0; wireCount < numWires; ++wireCount) {
                        unsigned inputCap = 0, inputSlew = 0;
                        unsigned tolerance = 0;
                        unsigned key = computeMinDelaySegment(charSegLength, load, outSlew,
                                                              inputCap, inputSlew, tolerance);
                        load = inputCap;
                        outSlew = inputSlew;
                        topology.addWireSegment(key); 
                }
        });
        
        topology.setInputSlew(outSlew);
        topology.setInputCap(load);
        
        //computeBranchingPoints(level, topology);
        
        _topologyForEachLevel[level-1] = topology;
}

unsigned HTreeBuilder::computeMinDelaySegment(unsigned length, unsigned targetLoad, unsigned targetOutSlew,
                                              unsigned &inputCap, unsigned &inputSlew, unsigned tolerance) const {
        unsigned minDelaySegKey = INVALID_KEY;
        unsigned minDelay = DELAY_MAX;
        
        unsigned numWireSegments = 0;
        
        unsigned minLoad = (targetLoad <= tolerance) ? (0) : (targetLoad - tolerance);
        unsigned maxLoad = targetLoad + tolerance;
        unsigned minOutSlew = (targetOutSlew <= tolerance) ? (0) : (targetOutSlew - tolerance);
        unsigned maxOutSlew = targetOutSlew + tolerance;

        bool forceBuffer = false;
        if (targetOutSlew == 1) {
                forceBuffer = true;                      
        }

        for (unsigned load = minLoad; load <= maxLoad; ++load) {
                for (unsigned outSlew = minOutSlew; outSlew <= maxOutSlew; ++outSlew) {
                        _techChar->forEachWireSegment(length, load, outSlew, 
                                                      [&] (unsigned key, const WireSegment& seg) {
                                if (seg.getOutputSlew() > _options->getMaxSlew()) {
                                       return;
                                }

                                //if (forceBuffer && !seg.isBuffered()) {
                                //        return;
                                //}

                                _techChar->reportSegment(key);
                                if (seg.getDelay() < minDelay) {
                                        minDelaySegKey = key;
                                        minDelay = seg.getDelay();
                                        inputCap = seg.getInputCap();
                                        inputSlew = seg.getInputSlew();
                                }
                                ++numWireSegments;
                        });
                }
        }
        
        if (numWireSegments == 0) {
                return computeMinDelaySegment(length, targetLoad, targetOutSlew, inputCap, inputSlew, tolerance + 1);
        }

        if (minDelaySegKey == INVALID_KEY) {
                std::cout << " [ERROR] Invalid key while computing min delay segment.\n";
                std::exit(1);        
        }
        
        std::cout << "    Tolerance: " << tolerance << " ";
        _techChar->reportSegment(minDelaySegKey);
        
        return minDelaySegKey;
}


void HTreeBuilder::computeBranchingPoints(unsigned level) {
        std::cout << " " << level << "... ";
        LevelTopology& topology = _topologyForEachLevel[level-1];

        if (level == 1) {
                Point<double> clockRoot(_sinkRegion.computeCenter());
                unsigned branchPtIdx1 = 
                        topology.addBranchingPoint(Point<double>(clockRoot.getX() - topology.getLength(), 
                                                                 clockRoot.getY()),
                                                   LevelTopology::NO_PARENT); 
                unsigned branchPtIdx2 = 
                        topology.addBranchingPoint(Point<double>(clockRoot.getX() + topology.getLength(), 
                                                                 clockRoot.getY()),
                                                   LevelTopology::NO_PARENT);
                std::vector<std::pair<float, float>> topLevelSinks;
                initTopLevelSinks(topLevelSinks);
                refineBranchingPointsWithClustering(topology, level, branchPtIdx1, branchPtIdx2, clockRoot,
                                                    topLevelSinks);
                return;
        }
        
        LevelTopology& parentTopology = _topologyForEachLevel[level-2];
        parentTopology.forEachBranchingPoint( [&] (unsigned idx, Point<double> clockRoot) {
                if (isHorizontal(level)) {
                        unsigned branchPtIdx1 = 
                                topology.addBranchingPoint(Point<double>(clockRoot.getX() - topology.getLength(), 
                                                                         clockRoot.getY()), idx);
                        unsigned branchPtIdx2 =  
                                topology.addBranchingPoint(Point<double>(clockRoot.getX() + topology.getLength(), 
                                                                         clockRoot.getY()), idx);
                        
                        std::vector<std::pair<float, float>> sinks;
                        computeBranchSinks(parentTopology, idx, sinks);
                        refineBranchingPointsWithClustering(topology, level, branchPtIdx1, branchPtIdx2, 
                                                            clockRoot, sinks);
                } else {
                        unsigned branchPtIdx1 = 
                                topology.addBranchingPoint(Point<double>(clockRoot.getX(), 
                                                                         clockRoot.getY() - topology.getLength()), idx); 
                        unsigned branchPtIdx2 =          
                                topology.addBranchingPoint(Point<double>(clockRoot.getX(), 
                                                                         clockRoot.getY() + topology.getLength()), idx);

                        std::vector<std::pair<float, float>> sinks;
                        computeBranchSinks(parentTopology, idx, sinks);
                        refineBranchingPointsWithClustering(topology, level, branchPtIdx1, branchPtIdx2, 
                                                            clockRoot, sinks);
                }
        });
}

void HTreeBuilder::initTopLevelSinks(std::vector<std::pair<float,float>>& sinkLocations) {
        sinkLocations.clear();
        _clock.forEachSink( [&] (const ClockInst& sink) {
                sinkLocations.emplace_back( (float) sink.getX() / _wireSegmentUnit, 
                                            (float) sink.getY() / _wireSegmentUnit );
        });
}

void HTreeBuilder::computeBranchSinks(LevelTopology& topology, unsigned branchIdx,
                                      std::vector<std::pair<float,float>>& sinkLocations) const {
        sinkLocations.clear();
        for (const Point<double>& point : topology.getBranchSinksLocations(branchIdx)) {
                sinkLocations.emplace_back(point.getX(), point.getY());
        }         
}

void HTreeBuilder::refineBranchingPointsWithClustering(LevelTopology& topology,
                                                       unsigned level,
                                                       unsigned branchPtIdx1,
                                                       unsigned branchPtIdx2, 
                                                       const Point<double>& rootLocation,
                                                       const std::vector<std::pair<float, float>>& sinks ) {
        //std::cout << "    Refining branch point locations\n";
        CKMeans::clustering clusteringEngine(sinks, rootLocation.getX(), rootLocation.getY());
        clusteringEngine.setPlotFileName("plot_" + std::to_string(level) + "_" + 
                                         std::to_string(branchPtIdx1) + "_" + 
                                         std::to_string(branchPtIdx2) + ".py");
        
        Point<double>& branchPt1 = topology.getBranchingPoint(branchPtIdx1);
        Point<double>& branchPt2 = topology.getBranchingPoint(branchPtIdx2);
        double targetDist = branchPt2.computeDist(rootLocation);
        //std::cout << "      R: "  << rootLocation << "\n"; 
        //std::cout << "      B1: " << branchPt1 << " B2: " << branchPt2 << "\n";
        //std::cout << "      D1: " << branchPt1.computeDist(rootLocation) 
        //          << " D2: " << branchPt2.computeDist(rootLocation) << "\n";
       
        std::vector<std::pair<float, float>> means;
        means.emplace_back(branchPt1.getX(), branchPt1.getY());
        //means.emplace_back(rootLocation.getX(), rootLocation.getY());
        means.emplace_back(branchPt2.getX(), branchPt2.getY());
        
        clusteringEngine.iterKmeans(1, means.size(), sinks.size()/means.size(), 0, means, 5);
        branchPt1 = Point<double>(means[0].first, means[0].second);
        branchPt2 = Point<double>(means[1].first, means[1].second);
        
        //std::cout << "      B1: " << branchPt1 << " B2: " << branchPt2 << "\n";
        //std::cout << "      D1: " << branchPt1.computeDist(rootLocation) 
        //          << " D2: " << branchPt2.computeDist(rootLocation) << "\n";
        
        std::vector<std::vector<unsigned>> clusters;
        clusteringEngine.getClusters(clusters);
        for (unsigned clusterIdx = 0; clusterIdx < clusters.size(); ++clusterIdx) {
                //std::cout << "    Cluster size: " << clusters[clusterIdx].size() << "\n";
                for (unsigned elementIdx = 0; elementIdx < clusters[clusterIdx].size(); ++elementIdx) {
                        unsigned sinkIdx = clusters[clusterIdx][elementIdx];
                        Point<double> sinkLoc(sinks[sinkIdx].first, sinks[sinkIdx].second);
                        if (clusterIdx == 0) {
                                topology.addSinkToBranch(branchPtIdx1, sinkLoc);
                        } else {
                                topology.addSinkToBranch(branchPtIdx2, sinkLoc);
                        }               
                }   
        }

        assert(std::abs(branchPt1.computeDist(rootLocation) - targetDist) < 0.001 && 
               std::abs(branchPt2.computeDist(rootLocation) - targetDist) < 0.001);
}


void HTreeBuilder::createClockSubNets() {
        std::cout << " Building clock sub nets...\n";
       
        DBU centerX = _sinkRegion.computeCenter().getX() * _wireSegmentUnit;
        DBU centerY = _sinkRegion.computeCenter().getY() * _wireSegmentUnit;
        ClockInst& rootBuffer = _clock.addClockBuffer("clkbuf_0", _options->getRootBuffer(), 
                                                             centerX, centerY); 
        Clock::SubNet& rootClockSubNet = _clock.addSubNet("clknet_0");
        rootClockSubNet.addInst(rootBuffer);

        // First level...
        LevelTopology &topLevelTopology = _topologyForEachLevel[0];
        topLevelTopology.forEachBranchingPoint( [&] (unsigned idx, Point<double> branchPoint) {
                SegmentBuilder builder("clkbuf_1_" + std::to_string(idx) + "_",
                                       "clknet_1_" + std::to_string(idx) + "_",
                                       _sinkRegion.computeCenter(), 
                                       branchPoint,
                                       topLevelTopology.getWireSegments(),
                                       _clock,
                                       rootClockSubNet,
                                       *_techChar,
                                       _wireSegmentUnit);
                builder.build();
                if (_topologyForEachLevel.size() == 1) {
                        builder.forceBufferInSegment(_options->getRootBuffer());
                        
                }
                topLevelTopology.setBranchDrivingSubNet(idx, *builder.getDrivingSubNet());                
        });
        
        // Others...
        for (unsigned levelIdx = 1; levelIdx < _topologyForEachLevel.size(); ++levelIdx) {
                LevelTopology& topology  = _topologyForEachLevel[levelIdx];
                topology.forEachBranchingPoint( [&] (unsigned idx, Point<double> branchPoint) {
                        unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
                        LevelTopology &parentTopology = _topologyForEachLevel[levelIdx - 1];
                        Point<double> parentPoint = parentTopology.getBranchingPoint(parentIdx);
                       
                        SegmentBuilder builder("clkbuf_" + std::to_string(levelIdx+1) + "_" + 
                                               std::to_string(idx) + "_",
                                               "clknet_" + std::to_string(levelIdx+1) + "_" + 
                                               std::to_string(idx) + "_",
                                               parentPoint, 
                                               branchPoint,
                                               topology.getWireSegments(),
                                               _clock,
                                               *parentTopology.getBranchDrivingSubNet(parentIdx),
                                               *_techChar,
                                               _wireSegmentUnit);
                        builder.build();
                        if (levelIdx == _topologyForEachLevel.size() - 1) {
                                builder.forceBufferInSegment(_options->getRootBuffer());
                        }
                        topology.setBranchDrivingSubNet(idx, *builder.getDrivingSubNet());                     
                });
        }

        // ---
        std::map<Point<double>, ClockInst*> mapLocationToSink;
        _clock.forEachSink( [&] (ClockInst& inst) {
                        Point<double> normLocation( (float) inst.getX() / _wireSegmentUnit,
                                                    (float) inst.getY() / _wireSegmentUnit);
                        mapLocationToSink[normLocation] = &inst;
                });
        
        LevelTopology& leafTopology = _topologyForEachLevel.back();
        unsigned levelIdx = _topologyForEachLevel.size() - 1;
        unsigned numSinks = 0;
        leafTopology.forEachBranchingPoint( [&] (unsigned idx, Point<double> branchPoint) {
                Clock::SubNet* subNet = leafTopology.getBranchDrivingSubNet(idx);
                subNet->setLeafLevel(true);
                
                const std::vector<Point<double>>& sinkLocs = leafTopology.getBranchSinksLocations(idx);
                for (const Point<double>& loc : sinkLocs) {
                        if (mapLocationToSink.find(loc) == mapLocationToSink.end()) {
                                std::cout << "Sink not found!\n";
                                std::exit(1);
                        }
                        
                        subNet->addInst(*mapLocationToSink[loc]);
                        ++numSinks;
                }  
        });

        std::cout << " Number of sinks covered: " << numSinks << "\n";
}

void HTreeBuilder::createSingleBufferClockNet() {
        std::cout << " Building single-buffer clock net...\n";
       
        DBU centerX = _sinkRegion.computeCenter().getX() * _wireSegmentUnit;
        DBU centerY = _sinkRegion.computeCenter().getY() * _wireSegmentUnit;
        ClockInst& rootBuffer = _clock.addClockBuffer("clkbuf_0", _options->getRootBuffer(), 
                                                             centerX, centerY); 
        Clock::SubNet& clockSubNet = _clock.addSubNet("clknet_0");
        clockSubNet.addInst(rootBuffer);

        _clock.forEachSink( [&] (ClockInst& inst) {
                clockSubNet.addInst(inst);                                               
        });
}

void HTreeBuilder::plotSolution() {
        std::ofstream file("plot.py");
        file << "import numpy as np\n"; 
        file << "import matplotlib.pyplot as plt\n";
        file << "import matplotlib.path as mpath\n";
        file << "import matplotlib.lines as mlines\n";
        file << "import matplotlib.patches as mpatches\n";
        file << "from matplotlib.collections import PatchCollection\n\n";

        _clock.forEachSink( [&] (const ClockInst& sink) {
                file << "plt.scatter(" << (double) sink.getX() / _wireSegmentUnit << ", " 
                     << (double) sink.getY() / _wireSegmentUnit << ")\n"; 
        });

        LevelTopology &topLevelTopology = _topologyForEachLevel.front();
        Point<double> topLevelBufferLoc = _sinkRegion.computeCenter();
        topLevelTopology.forEachBranchingPoint( [&] (unsigned idx, Point<double> branchPoint) {
                if (topLevelBufferLoc.getX() < branchPoint.getX()) {
                        file << "plt.plot([" 
                             << topLevelBufferLoc.getX() << ", " 
                             << branchPoint.getX() << "], ["
                             << topLevelBufferLoc.getY() << ", " 
                             << branchPoint.getY() << "], c = 'r')\n";
                } else {
                        file << "plt.plot([" 
                             << branchPoint.getX() << ", " 
                             << topLevelBufferLoc.getX() << "], ["
                             << branchPoint.getY() << ", " 
                             << topLevelBufferLoc.getY() << "], c = 'r')\n";
                }
        });        

        for (unsigned levelIdx = 1; levelIdx < _topologyForEachLevel.size(); ++levelIdx) {
                const LevelTopology& topology  = _topologyForEachLevel[levelIdx];
                topology.forEachBranchingPoint( [&] (unsigned idx, Point<double> branchPoint) {
                        unsigned parentIdx = topology.getBranchingPointParentIdx(idx);
                        Point<double> parentPoint = _topologyForEachLevel[levelIdx - 1].getBranchingPoint(parentIdx);
                        std::string color = "orange";
                        if (levelIdx % 2 == 0) {
                                color = "red";
                        } 

                        if (parentPoint.getX() < branchPoint.getX()) {
                                file << "plt.plot([" 
                                     << parentPoint.getX() << ", " 
                                     << branchPoint.getX() << "], ["
                                     << parentPoint.getY() << ", " 
                                     << branchPoint.getY() << "], c = '" << color <<  "')\n";
                        } else {
                                file << "plt.plot([" 
                                     << branchPoint.getX() << ", " 
                                     << parentPoint.getX() << "], ["
                                     << branchPoint.getY() << ", " 
                                     << parentPoint.getY() << "], c = '" << color << "')\n";
                        }
                                                
                });
        }

        file << "plt.show()\n";
        file.close();
}

void SegmentBuilder::build() {
        if (_root.getX() == _target.getX()) {
                buildVerticalConnection();
        } else if (_root.getY() == _target.getY()) { 
                buildHorizontalConnection();
        } else {
                buildLShapeConnection();
        }
}

void SegmentBuilder::buildVerticalConnection() {
        double length = std::abs(_root.getY() - _target.getY());
        bool isLowToHi = _root.getY() < _target.getY();
        
        double x = _root.getX();
        
        double currLength = 0;        
        for (unsigned wire = 0; wire < _techCharWires.size(); ++wire) {
                unsigned techCharWireIdx = _techCharWires[wire];
                const WireSegment& wireSegment = _techChar->getWireSegment(techCharWireIdx);
                unsigned wireSegLen = wireSegment.getLength();
                for (unsigned buffer = 0; buffer < wireSegment.getNumBuffers(); ++buffer) {
                        double location = wireSegment.getBufferLocation(buffer) * wireSegLen;
                        currLength += location;
                        double y = (isLowToHi) ? (_root.getY() + currLength) : (_root.getY() - currLength);
                        //std::cout << "B " << Point<double>(x, y) << " " 
                        //          << wireSegment.getBufferMaster(buffer) << "\n";
                        
                        _clock->addClockBuffer(_instPrefix + std::to_string(_numBuffers),
                                                     wireSegment.getBufferMaster(buffer),
                                                     x * _techCharDistUnit,
                                                     y * _techCharDistUnit);
                        ++_numBuffers;
                }
        }
        //std::cout << "currLength: " << currLength << " length: " << length << "\n";
}

void SegmentBuilder::buildHorizontalConnection() {
        double length = std::abs(_root.getX() - _target.getX());
        bool isLowToHi = _root.getX() < _target.getX();
        
        double y = _root.getY();
        
        double currLength = 0;        
        for (unsigned wire = 0; wire < _techCharWires.size(); ++wire) {
                unsigned techCharWireIdx = _techCharWires[wire];
                const WireSegment& wireSegment = _techChar->getWireSegment(techCharWireIdx);
                unsigned wireSegLen = wireSegment.getLength();
                for (unsigned buffer = 0; buffer < wireSegment.getNumBuffers(); ++buffer) {
                        double location = wireSegment.getBufferLocation(buffer) * wireSegLen;
                        currLength += location;
                        double x = (isLowToHi) ? (_root.getX() + currLength) : (_root.getX() - currLength);
                        //std::cout << "B " << Point<double>(x, y) << " " 
                        //          << wireSegment.getBufferMaster(buffer) << "\n";

                        _clock->addClockBuffer(_instPrefix + std::to_string(_numBuffers),
                                                  wireSegment.getBufferMaster(buffer),
                                                  x * _techCharDistUnit,
                                                  y * _techCharDistUnit);
                        ++_numBuffers;
                }
        }
        //std::cout << "currLength: " << currLength << " length: " << length << "\n";
}

void SegmentBuilder::buildLShapeConnection() {
        //std::cout << "L shape connection. Exiting...\n";
        //std::exit(1);
        double lengthX = std::abs(_root.getX() - _target.getX());
        double lengthY = std::abs(_root.getY() - _target.getY());
        bool isLowToHiX = _root.getX() < _target.getX();
        bool isLowToHiY = _root.getY() < _target.getY();

        //std::cout << "      Root: " << _root << " target: " << _target << "\n";
        double currLength = 0.0;
        for (unsigned wire = 0; wire < _techCharWires.size(); ++wire) {
                unsigned techCharWireIdx = _techCharWires[wire];
                const WireSegment& wireSegment = _techChar->getWireSegment(techCharWireIdx);
                unsigned wireSegLen = wireSegment.getLength();
                for (unsigned buffer = 0; buffer < wireSegment.getNumBuffers(); ++buffer) {
                        double location = wireSegment.getBufferLocation(buffer) * wireSegLen;
                        currLength += location;
                        
                        double x = std::numeric_limits<double>::max();
                        double y = std::numeric_limits<double>::max();
                        if (currLength < lengthX) {
                                y = _root.getY();
                                x = (isLowToHiX) ? (_root.getX() + currLength) : 
                                                   (_root.getX() - currLength);
                        } else {
                                x = _target.getX();
                                y = (isLowToHiY) ? (_root.getY() + (currLength - lengthX)) : 
                                                   (_root.getY() - (currLength - lengthX));
                        }

                        ClockInst& newBuffer = 
                                _clock->addClockBuffer(_instPrefix + std::to_string(_numBuffers),
                                                          wireSegment.getBufferMaster(buffer),
                                                          x * _techCharDistUnit,
                                                          y * _techCharDistUnit);
                        _drivingSubNet->addInst(newBuffer);
                        _drivingSubNet = &_clock->addSubNet(_netPrefix + std::to_string(_numBuffers));
                        _drivingSubNet->addInst(newBuffer);
                        
                        //std::cout << "      x: " << x << " y: " << y << "\n";
                        ++_numBuffers;
                }
        }
}

void SegmentBuilder::forceBufferInSegment(std::string master) {
        if (_numBuffers != 0) {
                return;
        }

        unsigned x = _target.getX();
        unsigned y = _target.getY();
        ClockInst& newBuffer = _clock->addClockBuffer(_instPrefix + "_f",
                                                             master,
                                                             x * _techCharDistUnit,
                                                             y * _techCharDistUnit);
        _drivingSubNet->addInst(newBuffer);
        _drivingSubNet = &_clock->addSubNet(_netPrefix + "_leaf");
        _drivingSubNet->addInst(newBuffer);
}

}
