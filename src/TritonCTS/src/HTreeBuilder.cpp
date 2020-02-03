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

void HTreeBuilder::run() {
        std::cout << " Generating H-Tree topology for net " << _clock.getName() << "...\n";
        std::cout << "    Tot. number of sinks: " << _clock.getNumSinks() << "\n";
       
        _clockTreeMaxDepth = _options->getClockTreeMaxDepth();
        _minInputCap = _techChar->getActualMinInputCap();
        _numMaxLeafSinks = _options->getNumMaxLeafSinks();
        _numMinLeafSinks = _options->getNumMinLeafSinks();
        _minLengthSinkRegion = _techChar->getMinSegmentLength() * 2;
        
        initSinkRegion();
        initTopology();
        computeCharWires();
        
        if (_topologyForEachLevel.size() < 1) {
                createSingleBufferClockNet();
                return;
        }
        
        computeBranchLocations();
        
        if (_options->getPlotSolution()) {
                plotSolution();
        }
        
        createClockSubNets();

        std::cout << " Clock topology of net \"" << _clock.getName() << "\" done.\n";
}

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

void HTreeBuilder::initTopology() {
        std::cout << " Initializing topology...\n";
        
        _branchingFactors = _options->getBranchingFactors();
        if (_branchingFactors.size() < 1) {
                computeBranchingFactors();
        } 

        _maxLevel = _branchingFactors.size();

        unsigned numSinks = _clock.getNumSinks();
        double height     = _sinkRegion.getHeight();
        double width      = _sinkRegion.getWidth();
        for (unsigned level = 1; level <= _maxLevel; ++level) {
                LevelTopology topology;
                
                unsigned fanout = _branchingFactors[level-1];
                topology.setBranchingFactor(fanout); 

                numSinks = std::ceil((double) numSinks / fanout);
                topology.setNumSinks(numSinks);       
                
                if (isVertical(level)) {
                        topology.setDirection(LevelTopology::VERTICAL);                
                        height /= fanout;
                } else {
                        topology.setDirection(LevelTopology::HORIZONTAL);                
                        width  /= fanout;
                }
                
                if (width < 2.0 || height < 2.0) {
                        std::cout << "    [WARNING] Sink region too small. Stopping at level "
                                  << level << "...\n";
                        break;
                } 

                topology.setHeight(height);
                topology.setWidth(width);

                unsigned length = computeNearestSegmentLength(width/2.0);
                if (isVertical(level)) {
                        length = computeNearestSegmentLength(height/2.0);
                }
                topology.createTopologyWire(length);

                for (unsigned wire = 1; wire < fanout/2; ++wire) {
                        unsigned length = computeNearestSegmentLength(width);
                        if (isVertical(level)) {
                                length = computeNearestSegmentLength(height);
                        }
                
                        topology.createTopologyWire(length);
                };

                _topologyForEachLevel.push_back(topology);
        }
        
        reportTopology();

        if (numSinks < _numMinLeafSinks) {
                std::cout << " Number of sinks on the sinks region is too small. Exiting...\n";
                std::exit(1);
        }

        if (numSinks >= _numMaxLeafSinks) {
                std::cout << " Number of sinks on the sinks region is too large. Exiting...\n";
                std::exit(1);
        }
}

void HTreeBuilder::computeBranchingFactors() {
        std::cout << " Using branching factor of 2...\n";
        const unsigned numSinks = _clock.getNumSinks();
        _maxLevel = 1;
        _branchingFactors.push_back(2);
        while (std::ceil((double) numSinks / pow(2, _maxLevel)) >= _numMaxLeafSinks) {
                _branchingFactors.push_back(2);
                ++_maxLevel; 
        }
}

void HTreeBuilder::reportTopology() const {
        std::cout << " -----------------------------------\n";
        std::cout << " " << std::setw(14) << "" << std::setw(7) << "Num " << "\n";
        std::cout << " "
                  << std::setw(7) << "Level" 
                  << std::setw(7) << "Fanout"
                  << std::setw(7) << "sinks" 
                  << std::setw(7) << "Width" 
                  << std::setw(7) << "Height" << "\n";
        std::cout << " -----------------------------------\n";
       
        std::cout << std::setprecision(4); 
        for (unsigned level = 1; level <= _maxLevel; ++level) {
                const LevelTopology& topology = _topologyForEachLevel[level-1];
                std::cout << " "
                          << std::setw(7) << level 
                          << std::setw(7) << topology.getBranchingFactor() 
                          << std::setw(7) << topology.getNumSinks() 
                          << std::setw(7) << topology.getWidth() 
                          << std::setw(7) << topology.getHeight() 
                          << "\n";
        }
        std::cout << " -----------------------------------\n";
}

inline 
unsigned HTreeBuilder::computeNearestSegmentLength(double length) {
        unsigned minLength = _minLengthSinkRegion / 2;

        if (minLength > 1 && length < minLength && _options->isFakeLutEntriesEnabled()) {
                unsigned minIndex = 1;
                _techChar->createFakeEntries(_minLengthSinkRegion, minIndex);
                _minLengthSinkRegion = 1;
        }
        
        unsigned segmentLength = std::round((double)length/(minLength))*minLength;
        return std::max<unsigned>(segmentLength, 1);
}

void HTreeBuilder::computeCharWires() {
        TopologyWire dummyRootWire(0);
        dummyRootWire.setOutputCap(_minInputCap);
        dummyRootWire.setOutputSlew(1);

        TopologyWire *prevWire = &dummyRootWire;
        for (unsigned level = 1; level <= _maxLevel; ++level) {
                std::cout << " Computing wires for level "<< level << "...\n";
                bool plotHeader = true;
                LevelTopology& topology = _topologyForEachLevel[level-1];
                topology.forEachTopologyWire( [&] (TopologyWire& wire) {
                        computeCharWire(wire, 
                                        prevWire->getOutputCap(), 
                                        prevWire->getOutputSlew());
                        prevWire = &wire;   
                        reportTopologyWire(wire, plotHeader);
                        plotHeader = false;
                });                                
        }
}

void HTreeBuilder::reportTopologyWire(TopologyWire& wire, bool plotHeader) const {
        if (plotHeader) {
                std::cout << " ---------------------------------------\n";
                std::cout << " "
                          << std::setw(6) << "Key" 
                          << std::setw(9) << "OutSlew"
                          << std::setw(6) << "Load"
                          << std::setw(9) << "Length"
                          << std::setw(9) << "isBuf" << "\n";
                std::cout << " ---------------------------------------\n";
        }
        for (unsigned wireKey: wire.allCharWires()) {
                const WireSegment& wire = _techChar->getWireSegment(wireKey);
                std::cout << std::setw(6) << wireKey 
                          << std::setw(9) << (unsigned) wire.getOutputSlew()
                          << std::setw(6) << (unsigned) wire.getLoad()
                          << std::setw(9) << (unsigned) wire.getLength()
                          << std::setw(9) << wire.isBuffered() << "\n";
        }
        std::cout << " ---------------------------------------\n";
};


void HTreeBuilder::computeCharWire(TopologyWire& wire, unsigned inputCap, unsigned inputSlew) {
        const unsigned SLEW_THRESHOLD = _options->getMaxSlew();
        const unsigned INIT_TOLERANCE = 1;

        unsigned segmentLength = wire.getLength();
        
        unsigned currLength = 0;
        for (unsigned charSegLength = _techChar->getMaxSegmentLength(); charSegLength >= 1; --charSegLength) {
                unsigned numWires = (segmentLength - currLength) / charSegLength;
                
                if (numWires < 1) {
                        continue;
                }

                currLength += numWires * charSegLength;
                for (unsigned wireCount = 0; wireCount < numWires; ++wireCount) {
                        unsigned outCap = 0, outSlew = 0;
                        unsigned key = computeMinDelaySegment(charSegLength, inputSlew, inputCap, 
                                                              SLEW_THRESHOLD, INIT_TOLERANCE, outSlew, outCap);
                        
                        inputCap = std::max(outCap, _minInputCap);
                        inputSlew = outSlew;
                        wire.addCharWire(key); 
                }

                if (currLength == segmentLength) {
                        break;
                }
        }
        
        wire.setOutputCap(inputCap);
        wire.setOutputSlew(inputSlew);
}

void HTreeBuilder::computeBranchLocations() {
        LevelTopology dummyRootTopology;
        dummyRootTopology.addBranchLocation(_sinkRegion.computeCenter());
        _clock.forEachSink( [&] (const ClockInst& sink) {
                dummyRootTopology.addSinkToBranch( 0, Point<double>((float) sink.getX() / _wireSegmentUnit, 
                                                                    (float) sink.getY() / _wireSegmentUnit));
        });

        LevelTopology* parentTopology = &dummyRootTopology;
        for (unsigned level = 1; level <= _maxLevel; ++level) {
                LevelTopology& topology = _topologyForEachLevel[level-1];
                topology.computeBranchLocations(parentTopology);
                parentTopology = &topology;
        }
}

void LevelTopology::computeBranchLocations(LevelTopology* parentTopology) {
        parentTopology->forEachBranchLocation([&] (unsigned idx, Point<double> loc) {
                Point<double> lowerLoc = loc;
                Point<double> upperLoc = loc;
                unsigned prevLowerBranch = NO_PARENT;
                unsigned prevUpperBranch = NO_PARENT;
                forEachTopologyWire([&] (unsigned wireIdx, TopologyWire& wire) {
                        if (isVertical()) {
                                lowerLoc.setY(lowerLoc.getY() - wire.getLength());                        
                                upperLoc.setY(upperLoc.getY() + wire.getLength());                        
                        } else {
                                lowerLoc.setX(lowerLoc.getX() - wire.getLength());                        
                                upperLoc.setX(upperLoc.getX() + wire.getLength());                        
                        }
                        unsigned lowerIdx = addBranchLocation(lowerLoc);       
                        addUpstreamBranchIdx(idx);
                        setParentBranchIdx(lowerIdx, prevLowerBranch);
                        prevLowerBranch = lowerIdx;
                        setTopologyWireIdx(lowerIdx, wireIdx);

                        unsigned upperIdx = addBranchLocation(upperLoc);       
                        addUpstreamBranchIdx(idx);
                        setParentBranchIdx(upperIdx, prevUpperBranch);
                        prevUpperBranch = upperIdx;
                        setTopologyWireIdx(upperIdx, wireIdx);
                });

                refineBranchLocations(idx, parentTopology->getBranchSinksLocations(idx), loc);
        });
}       

void LevelTopology::refineBranchLocations(unsigned parentIdx,
                                          const std::vector<Point<double>>& sinkLocs,
                                          const Point<double> rootLoc) {
        std::vector<std::pair<float,float>> points;
        for (unsigned sink = 0; sink < sinkLocs.size(); ++sink) {
                points.emplace_back(sinkLocs[sink].getX(), sinkLocs[sink].getY());
        }
        
        CKMeans::clustering ckmeans(points, rootLoc.getX(), rootLoc.getY());

        std::vector<std::pair<float, float>> means;
        std::map<unsigned, unsigned> meanToBranchIdx;
        forEachBranchLocation([&] (unsigned idx, Point<double> loc) {
                if (getUpstreamBranchIdx(idx) == parentIdx) {
                        means.emplace_back(loc.getX(), loc.getY());
                        meanToBranchIdx[means.size()-1] = idx;
                }
        });
        
        ckmeans.iterKmeans(1, means.size(), points.size()/means.size(), 0, means, 5);

        std::vector<std::vector<unsigned>> clusters;
        ckmeans.getClusters(clusters);

        for (unsigned clusterIdx = 0; clusterIdx < clusters.size(); ++clusterIdx) {
                unsigned branchIdx = meanToBranchIdx.at(clusterIdx);
                Point<double> newLocation(means[clusterIdx].first, means[clusterIdx].second);
                setBranchLocation(branchIdx, newLocation);

                for (unsigned elementIdx = 0; elementIdx < clusters[clusterIdx].size(); ++elementIdx) {
                        unsigned sinkIdx = clusters[clusterIdx][elementIdx];
                        Point<double> sinkLoc(points[sinkIdx].first, points[sinkIdx].second);
                        addSinkToBranch(branchIdx, sinkLoc);
                }   
        }
}

unsigned HTreeBuilder::computeMinDelaySegment(unsigned length, unsigned inputSlew, unsigned inputCap, unsigned slewThreshold, 
                                              unsigned tolerance, unsigned &outputSlew, unsigned &outputCap ) const {
        unsigned minKey      = std::numeric_limits<unsigned>::max();
        unsigned minDelay    = std::numeric_limits<unsigned>::max();
        unsigned minBufKey   = std::numeric_limits<unsigned>::max();
        unsigned minBufDelay = std::numeric_limits<unsigned>::max();
      
        for (unsigned load = 1; load <= _techChar->getMaxCapacitance(); ++load) {
                for (unsigned outSlew = 1; outSlew <= _techChar->getMaxSlew(); ++outSlew) {
                        _techChar->forEachWireSegment(length, load, outSlew,
                                [&] (unsigned key, const WireSegment& seg) {
                                        //if (seg.getInputCap() != inputCap || seg.getInputSlew() != inputSlew) {
                                        //        return;
                                        //}

                                        //std::cout << "[" << length << ", " << load << ", " << outSlew << "] =  "
                                        //          << (unsigned) seg.getLength() << ", " << (unsigned) seg.getLoad() << ", " << (unsigned) seg.getOutputSlew()
                                        //          << ", " << key << ", " << _techChar->computeKey(length, load, outSlew) << "\n";                                       
                                        
                                        if ( std::abs( (int) seg.getInputCap() - (int) inputCap ) > tolerance || 
                                             std::abs( (int) seg.getInputSlew() - (int) inputSlew ) > tolerance ) {
                                                return;
                                        }

                                        //std::cout << "[" << length << ", " << load << ", " << outSlew << "] =  "
                                        //          << (unsigned) seg.getLength() << ", " << (unsigned) seg.getLoad() << ", " << (unsigned) seg.getOutputSlew()
                                        //          << ", " << key << ", " << _techChar->computeKey(length, load, outSlew) << "\n";                                       
                                        

                                        if (seg.getDelay() < minDelay) {
                                                minDelay = seg.getDelay();
                                                minKey = key;
                                        }

                                        if (seg.isBuffered() && seg.getDelay() < minBufDelay) {
                                                minBufDelay = seg.getDelay();
                                                minBufKey = key;
                                        }
                                });                
                }
        }

        const unsigned MAX_TOLERANCE = 10;
        if (inputSlew >= slewThreshold) {
                if (minBufKey < std::numeric_limits<unsigned>::max()) {
                        const WireSegment& bestBufSegment = _techChar->getWireSegment(minBufKey);
                        outputSlew = bestBufSegment.getOutputSlew();
                        outputCap  = bestBufSegment.getLoad();      
                        return minBufKey;
                } else if (tolerance < MAX_TOLERANCE) {
                        //std::cout << "    Could not find a buffered segment for [" << length << ", " 
                        //          << inputSlew << ", " << inputCap << "]... ";
                        //std::cout << "Increasing tolerance\n";
                        return computeMinDelaySegment(length, inputSlew, inputCap, slewThreshold, 
                                                      tolerance + 1, outputSlew, outputCap );
                }
        }

        if (minKey == std::numeric_limits<unsigned>::max()) {
                //std::cout << "    Could not find segment for [" << length << ", " 
                //          << inputSlew << ", " << inputCap << "]... ";
                //std::cout << "Increasing tolerance\n";
                return computeMinDelaySegment(length, inputSlew, inputCap, slewThreshold, 
                                              tolerance + 1, outputSlew, outputCap );
        }

        const WireSegment& bestSegment = _techChar->getWireSegment(minKey);
        outputSlew = std::max( (unsigned) bestSegment.getOutputSlew(), inputSlew + 1 );
        outputCap  = bestSegment.getLoad();

        return minKey;
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
        topLevelTopology.forEachBranchLocation( [&] (unsigned idx, Point<double> branchPoint) {
                Clock::SubNet* drivingSubNet = &rootClockSubNet;
                Point<double> source = _sinkRegion.computeCenter();
                unsigned parentBranchIdx = topLevelTopology.getParentBranchIdx(idx);
                if (parentBranchIdx != LevelTopology::NO_PARENT) {
                        drivingSubNet = topLevelTopology.getBranchDrivingSubNet(parentBranchIdx);
                        source = topLevelTopology.getBranchLocation(parentBranchIdx);
                }                

                SegmentBuilder builder("clkbuf_1_" + std::to_string(idx) + "_",
                                       "clknet_1_" + std::to_string(idx) + "_",
                                       source, 
                                       branchPoint,
                                       topLevelTopology.getTopologyWire(idx).allCharWires(),
                                       _clock,
                                       *drivingSubNet,
                                       *_techChar,
                                       _wireSegmentUnit);
                builder.build();
                if (_topologyForEachLevel.size() == 1) {
                        builder.forceBufferInSegment(_options->getRootBuffer());
                        
                }
                topLevelTopology.setBranchDrivingSubNet(idx, *builder.getDrivingSubNet());                
        });
        //return; 
        // Others...
        for (unsigned levelIdx = 1; levelIdx < _topologyForEachLevel.size(); ++levelIdx) {
                LevelTopology& topology  = _topologyForEachLevel[levelIdx];
                topology.forEachBranchLocation( [&] (unsigned idx, Point<double> branchPoint) {
                        LevelTopology &parentTopology = _topologyForEachLevel[levelIdx - 1];
                        
                        Clock::SubNet* drivingSubNet = nullptr;
                        Point<double> source(0, 0);
                        unsigned parentBranchIdx = topology.getParentBranchIdx(idx);
                        if (parentBranchIdx != LevelTopology::NO_PARENT) {
                                source = topology.getBranchLocation(parentBranchIdx);
                                drivingSubNet = topology.getBranchDrivingSubNet(parentBranchIdx);
                        } else {
                                parentBranchIdx = topology.getUpstreamBranchIdx(idx);
                                source = parentTopology.getBranchLocation(parentBranchIdx);
                                drivingSubNet = parentTopology.getBranchDrivingSubNet(parentBranchIdx);
                        }

                        SegmentBuilder builder("clkbuf_" + std::to_string(levelIdx+1) + "_" + 
                                               std::to_string(idx) + "_",
                                               "clknet_" + std::to_string(levelIdx+1) + "_" + 
                                               std::to_string(idx) + "_",
                                               source, 
                                               branchPoint,
                                               topology.getTopologyWire(idx).allCharWires(),
                                               _clock,
                                               *drivingSubNet,
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
        leafTopology.forEachBranchLocation( [&] (unsigned idx, Point<double> branchPoint) {
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
        topLevelTopology.forEachBranchLocation( [&] (unsigned idx, Point<double> branchPoint) {
                Point<double> parentLoc = topLevelBufferLoc;
                unsigned parentIdx = topLevelTopology.getParentBranchIdx(idx);
                if (parentIdx != LevelTopology::NO_PARENT) {
                        parentLoc = topLevelTopology.getBranchLocation(parentIdx);
                }                
                        
                if (parentLoc.getX() < branchPoint.getX()) {
                        file << "plt.plot([" 
                             << parentLoc.getX() << ", " 
                             << branchPoint.getX() << "], ["
                             << parentLoc.getY() << ", " 
                             << branchPoint.getY() << "], c = 'r')\n";
                } else {
                        file << "plt.plot([" 
                             << branchPoint.getX() << ", " 
                             << parentLoc.getX() << "], ["
                             << branchPoint.getY() << ", " 
                             << parentLoc.getY() << "], c = 'r')\n";
                }
        });        
       
        for (unsigned levelIdx = 1; levelIdx < _topologyForEachLevel.size(); ++levelIdx) {
                const LevelTopology& topology  = _topologyForEachLevel[levelIdx];
                topology.forEachBranchLocation( [&] (unsigned idx, Point<double> branchPoint) {
                        unsigned parentIdx = topology.getParentBranchIdx(idx);
                        Point<double> parentPoint(0, 0);
                        if (parentIdx != LevelTopology::NO_PARENT) {
                                parentPoint = topology.getBranchLocation(parentIdx);
                        } else {
                                unsigned upstreamBranchIdx = topology.getUpstreamBranchIdx(idx);
                                parentPoint = _topologyForEachLevel[levelIdx - 1].getBranchLocation(upstreamBranchIdx);
                        }

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
