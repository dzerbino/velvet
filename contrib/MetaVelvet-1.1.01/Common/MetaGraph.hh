#ifndef _META_GRAPH_HH_
#define _META_GRAPH_HH_
#include <math.h>
#include "../VelvetAPI/VelvetGraph.hh"
#include "InfiniteLoopChecker.hh"

#define META_GRAPH_MAX_NUM_COVERAGE_PEAKS 100
#define META_GRAPH_MAX_NUM_INTER_LOOPS 100
#define META_GRAPH_MASK_UNVISITED 0
#define META_GRAPH_MASK_NOW_VISITING 1
#define META_GRAPH_MASK_PRIMARY 2
#define META_GRAPH_MASK_SECONDARY -1
#define META_GRAPH_MASK_DELETED -2
#define META_GRAPH_LN2 0.693147
#define META_GRAPH_LONG_NODE_CUTOFF 50
#define META_GRAPH_UNIQUE_PROBABILITY_CUTOFF 5

using namespace std;

class MetaGraph : public VelvetGraph {
  int numCoveragePeaks;
  double expectedCoverages[META_GRAPH_MAX_NUM_COVERAGE_PEAKS];
  int* subgraphMask;
  void scaffoldingWithMultiPeakMode( boolean* flagMatePair, bool flagScaffolding, double maxChimeraRate, double repeatCoverageSD, bool flagDiscardChimera );
  void eliminateNullNodes();
  void extractNextSubgraph( IDnum nodeID );
  IDnum getUnvisitedNodeID() const;
  void dfs( IDnum currentIndex );
  bool isChimeraSubgraph( double primaryCoverage, double secondaryCoverage, double maxChimeraRate ) const;
  bool isPrimarySubgraph( double primaryCoverage, double secondaryCoverage ) const;
  int splitRepeats( double primaryCoverage, double secondaryCoverage, double repeatCoverageSD );
  bool isLocalRepeatStructure( Node* node ) const;
  bool isRepeatCoverageCondition( Node* node, double repeatCoverageSD ) const;
  double getNearestPeak( double coverage ) const ;
  bool pushNeighboursInterRepeat( Node* inNode, Node* repNode, Node* outNode );
  void absorbExtensionInterRepeat( Node * node, Node * extension );
  void adjustShortReadsInterRepeat( Node * target, Node * source );
  void forceSeparateChimericSubgraph( double primaryCoverage, double secondaryCoverage );
  void identifyUniqueNodesSubgraph( double primaryCoverage );
  bool isUniqueSubgraph(Node * node, double expCovSubgraph);
  void readCoherentSubgraph( double primaryCoverage );
  bool trimLongReadTipsSubgraph();
  void changeSubgraphMask( int fromID, int toID );
  long getSubgraphMaskSize() const;
  bool checkLongReadExistence() const;
  void resetNodeFlags();
  void resetUniqueness();
  IDnum getNumActiveNodes() const;
public:
  MetaGraph( const string& seqFileName, const string& roadmapFileName, const string& pregraphFileName, bool flagReadTracking, int accelerationBits ) 
    : VelvetGraph( seqFileName, roadmapFileName, pregraphFileName, flagReadTracking, accelerationBits ){}
  MetaGraph( const string& seqFileName, const string& graphFileName ) : VelvetGraph( seqFileName, graphFileName ){}
  void setExpectedCoverages( int num, const double* covs );
  void scaffolding( boolean* flagMatePair, bool flagScaffolding, double maxChimeraRate, double repeatCoverageSD, bool flagDiscardChimera );
  void showExpectedCoverages() const;
};

#endif // _META_GRAPH_HH_
