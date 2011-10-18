#include "MetaGraph.hh"

void MetaGraph::setExpectedCoverages( int num, const double* covs ){
  numCoveragePeaks = num;
  for( int i=0 ; i<META_GRAPH_MAX_NUM_COVERAGE_PEAKS ; ++i ){ 
    expectedCoverages[i] = covs[i]; 
  }
}

void MetaGraph::showExpectedCoverages() const {
  for( int i=0 ; i<numCoveragePeaks ; ++i ){
    cout << "[MetaGraph] " << (i+1) << "-th coverage peak = " << expectedCoverages[i] << endl;
  }
}


void MetaGraph::scaffolding( boolean* flagMatePair, bool flagScaffolding, double maxChimeraRate, double repeatCoverageSD, bool flagDiscardChimera ){
  if( numCoveragePeaks == 1 ){
    cout << "[MetaGraph] " << " === Scaffolding with single peak mode ===" << endl;
    return scaffoldingWithSinglePeakMode( expectedCoverages[0], flagMatePair, flagScaffolding );
  }
  cout << "[MetaGraph] " << " === Scaffolding with multi-peak mode ===" << endl;
  scaffoldingWithMultiPeakMode( flagMatePair, flagScaffolding, maxChimeraRate, repeatCoverageSD, flagDiscardChimera );
}

void MetaGraph::scaffoldingWithMultiPeakMode( boolean* flagMatePair, bool flagScaffolding, double maxChimeraRate, double repeatCoverageSD, bool flagDiscardChimera ){
  cout << "[MetaGraph] " << " === Create read paring array ===" << endl;
  createReadPairingArray(sequences);

  cout << "[MetaGraph] " << " === Detach dubious reads ===" << endl;
  detachDubiousReads( sequences, dubious );

  cout << "[MetaGraph] " << " === Activate gap markers ===" << endl;
  activateGapMarkers(graph);

  cout << "[MetaGraph] " << " === Split into subgraphs ===" << endl;
  bool flagLongRead = checkLongReadExistence();
  subgraphMask = callocOrExit( getSubgraphMaskSize() , int );
  eliminateNullNodes();

  int roundID = 0;
  int primaryPeakID = 0;
  bool flagVisitedAllNodes = false;
  while( !flagVisitedAllNodes ){
    double primaryCoverage = expectedCoverages[primaryPeakID];
    double secondaryCoverage = expectedCoverages[primaryPeakID+1];
    cout << "[MetaGraph] " << (roundID+1) << "-th round" << endl
	 << "[MetaGraph] " << "Primary coverage   = " << primaryCoverage << endl
	 << "[MetaGraph] " << "Secondary coverage = " << secondaryCoverage << endl;
    if( secondaryCoverage > 0 ){
      ++primaryPeakID;
    }
    ++roundID;

    while( IDnum nodeID = getUnvisitedNodeID() ){
      extractNextSubgraph( nodeID );
      if( isChimeraSubgraph( primaryCoverage, secondaryCoverage, maxChimeraRate ) ){
	splitRepeats( primaryCoverage, secondaryCoverage, repeatCoverageSD );
	eliminateNullNodes();
	if( flagDiscardChimera ){
	  changeSubgraphMask( META_GRAPH_MASK_NOW_VISITING, META_GRAPH_MASK_DELETED );
	} else {
	  forceSeparateChimericSubgraph( primaryCoverage, secondaryCoverage );
	}
      } else if( isPrimarySubgraph( primaryCoverage, secondaryCoverage ) ){
	changeSubgraphMask( META_GRAPH_MASK_NOW_VISITING, META_GRAPH_MASK_PRIMARY );
      }
    }

    resetNodeFlags();
    cout << "[MetaGraph] " << " === Subgraph contiging ===" << endl;
    identifyUniqueNodesSubgraph( primaryCoverage ); 
    
    if( flagLongRead ){
      cout << "[MetaGraph] " << " === Subgraph Rock Band ===" << endl;
      readCoherentSubgraph( primaryCoverage );
    }
    
    cout << "[MetaGraph] " << " === Subgraph Scaffolding ===" << endl;
    for ( int pebbleRounds=pairedCategories(sequences)+1 ; pebbleRounds>0 ; --pebbleRounds){
      exploitShortReadPairs( graph, sequences, dubious, flagMatePair, flagScaffolding );
    }
    cout << "[MetaGraph] " << " === Done: Subgraph contiging & scaffolding === " << endl << endl;
    
    eliminateNullNodes();
    resetUniqueness();
    changeSubgraphMask( META_GRAPH_MASK_PRIMARY, META_GRAPH_MASK_DELETED );
    changeSubgraphMask( META_GRAPH_MASK_SECONDARY, META_GRAPH_MASK_UNVISITED );
    if( getUnvisitedNodeID() == 0 ){
      flagVisitedAllNodes = true;
    }
  }
  free(subgraphMask);	
}



void MetaGraph::eliminateNullNodes(){
  for( IDnum i = 0; i<getSubgraphMaskSize() ; ++i ) {
    Node* node = getNodeInGraph( graph, i-nodeCount(graph) );
    if( node == NULL || getNodeID(node)==0 ){
      subgraphMask[i] = META_GRAPH_MASK_DELETED;
    }
  }
}

void MetaGraph::extractNextSubgraph( IDnum nodeID ){
  dfs( nodeID );
  dfs( -1 * nodeID );
}

IDnum MetaGraph::getUnvisitedNodeID() const {
  long maskSize = getSubgraphMaskSize();
  for( long i=0 ; i<maskSize ; ++i ){
    if( i == (maskSize-1)/2 ){
      continue;
    }
    if( subgraphMask[i] == META_GRAPH_MASK_UNVISITED ){
      return i - (maskSize-1) / 2;
    }
  }
  return NULL_IDX;
}

void MetaGraph::dfs( IDnum currentIndex ){
  if( subgraphMask[currentIndex + nodeCount(graph)] == META_GRAPH_MASK_UNVISITED ){
    subgraphMask[currentIndex + nodeCount(graph)] = META_GRAPH_MASK_NOW_VISITING; 
    Arc *activeArc = NULL;
    for( activeArc=getArc(getNodeInGraph(graph, currentIndex)) ; activeArc!=NULL ; activeArc=getNextArc(activeArc) ){
      long nextIndex = getNodeID(getDestination(activeArc));
      if( subgraphMask[nextIndex] == META_GRAPH_MASK_UNVISITED ){
	dfs( nextIndex );
	dfs( -1 * nextIndex );
      }
    }
  }
}

bool MetaGraph::isChimeraSubgraph( double primaryCoverage, double secondaryCoverage, double maxChimeraRate ) const {
  double totalPrimaryLength = 0.0;
  double totalSecondaryLength = 0.0;
  for ( long i=1 ; i<=nodeCount(graph) ; ++i ){
    if( subgraphMask[i + nodeCount(graph)] == META_GRAPH_MASK_NOW_VISITING ){
      Node* node = getNodeInGraph( graph, i );
      if( node == NULL )
	continue;
      double curCoverage = VUtils::getNodeDensity( node );
      if ( fabs(curCoverage-primaryCoverage) <= fabs(curCoverage-secondaryCoverage) ){
	totalPrimaryLength += getNodeLength( node );
      } else {
	totalSecondaryLength += getNodeLength( node );
      }
    }
  }
  double primaryFraction = totalPrimaryLength / (totalPrimaryLength + totalSecondaryLength);
  double secondaryFraction = totalSecondaryLength / (totalPrimaryLength + totalSecondaryLength);
  if( primaryFraction>maxChimeraRate and secondaryFraction>maxChimeraRate ){
    return true;
  }
  return false;
}

bool MetaGraph::isPrimarySubgraph( double primaryCoverage, double secondaryCoverage ) const {
  double totalPrimaryLength = 0.0;
  double totalSecondaryLength = 0.0;
  for ( long i=1 ; i<=nodeCount(graph) ; ++i ){
    if( subgraphMask[i + nodeCount(graph)] == META_GRAPH_MASK_NOW_VISITING ){
      Node* node = getNodeInGraph( graph, i );
      if( node == NULL )
	continue;
      double curCoverage = VUtils::getNodeDensity( node );
      if( fabs(curCoverage - primaryCoverage) <= fabs(curCoverage-secondaryCoverage) ){
	totalPrimaryLength += getNodeLength( node );
      } else {
	totalSecondaryLength += getNodeLength( node );
      }
    }
  }
  double primaryFraction = totalPrimaryLength / (totalPrimaryLength + totalSecondaryLength);
  double secondaryFraction = totalSecondaryLength / (totalPrimaryLength + totalSecondaryLength);
  if( primaryFraction > secondaryFraction ){
    return true;
  }
  return false;
}

int MetaGraph::splitRepeats( double primaryCoverage, double secondaryCoverage, double repeatCoverageSD ){
  int numInterRepeats = 0;
  vector<Node*> inNodes, outNodes;
  resetNodeFlags();
  for ( IDnum i=0; i<nodeCount(graph) ; ++i ){
    Node* node = getNodeInGraph( graph, i+1 );
    if( getNodeID(node) == NULL_IDX ){
      continue;
    }
    if( subgraphMask[getNodeID(node)+nodeCount(graph)] != META_GRAPH_MASK_NOW_VISITING ){
      continue;
    }

    if( isLocalRepeatStructure(node) and isRepeatCoverageCondition(node, repeatCoverageSD) ){
      ++numInterRepeats;
      vector<Node*> inNodes  = VUtils::getInNodes( node );
      vector<Node*> outNodes = VUtils::getOutNodes( node );
      Node* primaryInNode    = VUtils::getMaxDensityNode( inNodes );
      Node* secondaryInNode  = VUtils::getMinDensityNode( inNodes );
      Node* primaryOutNode   = VUtils::getMaxDensityNode( outNodes );
      Node* secondaryOutNode = VUtils::getMinDensityNode( outNodes );
      cout << endl
	   << "[MetaGraph] REPEAT Node:        ID = " << getNodeID(node)
	   << "\tCov = " << VUtils::getNodeDensity(node)  << "\tLength = " << getNodeLength(node) << endl
	   << "[MetaGraph] Primary IN Node:    ID = " << getNodeID(primaryInNode)  
	   << "\tCov = " << VUtils::getNodeDensity(primaryInNode)  << "\tLength = " << getNodeLength(primaryInNode) << endl
	   << "[MetaGraph] Primary OUT Node:   ID = " << getNodeID(primaryOutNode)
	   << "\tCov = " << VUtils::getNodeDensity(primaryOutNode)  << "\tLength = " << getNodeLength(primaryOutNode) << endl
	   << "[MetaGraph] Secondary IN Node:  ID = " << getNodeID(secondaryInNode)
	   << "\tCov = " << VUtils::getNodeDensity(secondaryInNode) << "\tLength = " << getNodeLength(secondaryInNode) << endl
	   << "[MetaGraph] Secondary OUT Node: ID = " << getNodeID(secondaryOutNode)
	   << "\tCov = " << VUtils::getNodeDensity(secondaryOutNode) << "\tLength = " << getNodeLength(secondaryOutNode) << endl;
      setNodeStatus(node, true); setUniqueness(node, false);
      if( getNodeID(primaryInNode) != getNodeID(primaryOutNode) ){ // Taking into account the knock-turn exception
	setNodeStatus(primaryInNode, true); setUniqueness(primaryInNode, true); setNodeStatus(primaryOutNode, true); setUniqueness(primaryOutNode, true);
	pushNeighboursInterRepeat( primaryInNode, node, primaryOutNode );
      } else {
	
      }
      if( getNodeID(secondaryInNode) != getNodeID(secondaryOutNode) ){ // Taking into account the knock-turn exception
	setNodeStatus(secondaryOutNode, true); setUniqueness(secondaryOutNode, true); setNodeStatus(secondaryOutNode, true); setUniqueness(secondaryOutNode, true);
	pushNeighboursInterRepeat( secondaryInNode, node, secondaryOutNode );
      }
    }
  }
  resetNodeFlags();
  return numInterRepeats;
}

bool MetaGraph::isLocalRepeatStructure( Node * node ) const {
  if( !(VUtils::getNumOutArcs(node) == 2) ){
    return false;
  }
  if( !(VUtils::getNumInArcs(node) == 2) ){
    return false;
  }
  return true;
}

bool MetaGraph::isRepeatCoverageCondition( Node* node, double repeatCoverageSD ) const {
  vector<Node*> inNodes  = VUtils::getInNodes( node );
  vector<Node*> outNodes = VUtils::getOutNodes( node );
  double inPeak1  = getNearestPeak( VUtils::getMaxNodeDensity( inNodes ) );
  double inPeak2  = getNearestPeak( VUtils::getMinNodeDensity( inNodes ) );
  double outPeak1 = getNearestPeak( VUtils::getMaxNodeDensity( outNodes ) );
  double outPeak2 = getNearestPeak( VUtils::getMinNodeDensity( outNodes ) );
  if( inPeak1==outPeak1 and inPeak2==outPeak2 and inPeak1>inPeak2 ){
    double aveDensity = ( VUtils::getTotalNodeDensity(inNodes) + VUtils::getTotalNodeDensity(outNodes) ) / 2.0;
    if( VUtils::getNodeDensity(node)<=aveDensity*(1.0+repeatCoverageSD) and VUtils::getNodeDensity(node)>=aveDensity*(1.0-repeatCoverageSD) ){
      return true;
    }
  }
  return false;
}

double MetaGraph::getNearestPeak( double coverage ) const {
  double nearestPeak = expectedCoverages[0];
  double nearestDiff = fabs( coverage - nearestPeak );
  for( int i=1 ; i<numCoveragePeaks ; ++i ){
    double curPeak = expectedCoverages[i];
    double curDiff = fabs( coverage - curPeak );
    if( curDiff < nearestDiff ){
      nearestPeak = curPeak;
      nearestDiff = curDiff;
    }
  }
  return nearestPeak;
}

bool MetaGraph::pushNeighboursInterRepeat( Node* inNode, Node* repNode, Node* outNode ){
  NodeList* path = allocateNodeList();
  path->node = repNode;
  path->next = allocateNodeList();
  path->next->node = outNode;
  path->next->next = NULL;
  while( path ){
    Node* candidate = path->node;
    NodeList* tmpPath = path->next;
    deallocateNodeList( path );
    path = tmpPath;
    if( getUniqueness(candidate) ){
      concatenateReadStarts( inNode, candidate, graph );
      concatenateLongReads( inNode, candidate, graph );
      absorbExtensionInterRepeat( inNode, candidate );
      for( int cat = 0 ; cat<CATEGORIES ; ++cat ){
	incrementVirtualCoverage( inNode, (Category)cat, getVirtualCoverage(candidate, (Category)cat) );
	incrementOriginalVirtualCoverage( inNode, (Category)cat, getOriginalVirtualCoverage(candidate, (Category)cat) );
      }
      destroyNode( candidate, graph );
      return true;
    } else {
      adjustShortReadsInterRepeat( inNode, candidate );
      adjustLongReads( inNode, getNodeLength(candidate) );
      absorbExtensionInterRepeat( inNode, candidate );
    }
  }
  return false;
}

/*
bool MetaGraph::pushNeighboursInterRepeat( Node * node, Node * nodeInterRepeat, Node * oppositeNode ){
  Node *candidate;
  Node *lastCandidate = NULL;
  Category cat;
  NodeList *path, *tmp;
	
  // Make path (= NodeList of node and oppositeNode)
  path = allocateNodeList();
  path->node = nodeInterRepeat;
  path->next = allocateNodeList();
  path->next->node = oppositeNode;
  path->next->next = NULL;
  // Original
  while (path) {
    candidate = path->node;
    tmp = path->next;
    deallocateNodeList(path);
    path = tmp;
    
    if (getUniqueness(candidate)) {
      concatenateReadStarts(node, candidate, graph);
      concatenateLongReads(node, candidate, graph);
      absorbExtensionInterRepeat(node, candidate);
      // Read coverage
      for (cat = 0; cat < CATEGORIES; cat++) {
	incrementVirtualCoverage(node, cat, getVirtualCoverage(candidate, cat));
	incrementOriginalVirtualCoverage(node, cat, getOriginalVirtualCoverage(candidate, cat));
      }

      printf("\tConcatenated InNode %d -- OutNode %d\n", getNodeID(node), getNodeID(candidate));
      destroyNode(candidate, graph);
      return true;
    } else {
      adjustShortReadsInterRepeat(node, candidate);
      adjustLongReads( node, getNodeLength(candidate) );
      absorbExtensionInterRepeat(node, candidate);
      lastCandidate = candidate;
    }
  }
  return false;
}
*/

void MetaGraph::absorbExtensionInterRepeat( Node * node, Node * extension ){
  appendNodeGaps( node, extension, graph );
  appendDescriptors( node, extension ); // cause error ... ?
  while (getArc(node) != NULL){
    destroyArc( getArc(node), graph );
  }
  // Create new
  for( Arc* arc = getArc(extension) ; arc != NULL ; arc = getNextArc(arc) ){
    createAnalogousArc( node, getDestination(arc), arc, graph );
  }
}

void MetaGraph::adjustShortReadsInterRepeat( Node * target, Node * source ){
  if( !readStartsAreActivated(graph) ){
    return;
  }
  ShortReadMarker* targetArray = getNodeReads( getTwinNode(target), graph );
  IDnum targetLength = getNodeReadCount( getTwinNode(target), graph );
  Coordinate nodeLength = getNodeLength( source );
  
  for( IDnum index=0 ; index<targetLength ; ++index ){
    ShortReadMarker* marker = getShortReadMarkerAtIndex( targetArray, index );
    Coordinate position = getShortReadMarkerPosition( marker );
    position += nodeLength;
    setShortReadMarkerPosition( marker, position );
  }
}

void MetaGraph::forceSeparateChimericSubgraph( double primaryCoverage, double secondaryCoverage ){
  for( IDnum i=0; i<getSubgraphMaskSize() ; ++i){
    if( subgraphMask[i] == META_GRAPH_MASK_NOW_VISITING ){
      Node* node = getNodeInGraph( graph, i - nodeCount(graph) );
      if (node == NULL)
	continue;
      double curCoverage = VUtils::getNodeDensity(node);
      if( fabs(curCoverage - primaryCoverage) <= fabs(curCoverage - secondaryCoverage)){
	subgraphMask[i] = META_GRAPH_MASK_PRIMARY;
      } else {
	subgraphMask[i] = META_GRAPH_MASK_SECONDARY;
      }
    }
  }
}

void MetaGraph::identifyUniqueNodesSubgraph( double primaryCoverage ){
  for( IDnum i=0; i<nodeCount(graph) ; ++i ){
    Node* node = getNodeInGraph( graph, i+1 );
    if( node != NULL ){
      if( subgraphMask[i+1+nodeCount(graph)] == META_GRAPH_MASK_PRIMARY ){
	setUniqueness( node, isUniqueSubgraph(node, primaryCoverage) );
      } else {
	setUniqueness(node, false);
      }
    }
  }
}

bool MetaGraph::isUniqueSubgraph( Node * node, double expCoverage ){
  int nodeLength = getNodeLength(node);
  if (nodeLength > META_GRAPH_LONG_NODE_CUTOFF) {
    double nodeDensity = VUtils::getNodeDensity( node );
    double probability = -1 * META_GRAPH_LN2/2 + nodeLength/(2*expCoverage) * (expCoverage*expCoverage - nodeDensity*nodeDensity/2);
    return probability > META_GRAPH_UNIQUE_PROBABILITY_CUTOFF;
  }
  return false;
}

void MetaGraph::readCoherentSubgraph( double primaryCoverage ){
  RecycleBin* listMemory = newRecycleBin(sizeof(PassageMarkerList), 100000);
  if( !trimLongReadTipsSubgraph() ){
    destroyRecycleBin(listMemory);
    return;
  }
  IDnum previousNodeCount = 0;
  int checkModified = -1;
  while( previousNodeCount != getNumActiveNodes() ){
    previousNodeCount = getNumActiveNodes();
    for( IDnum nodeIndex=1 ; nodeIndex<=nodeCount(graph) ; nodeIndex++ ){
      Node* node = getNodeInGraph(graph, nodeIndex);
      if( node == NULL || !getUniqueness(node) ){
	continue;
      }
      while( uniqueNodesConnect(node) ){
	node = bypass();
      }
      node = getTwinNode(node);
      while( uniqueNodesConnect(node) ){
	node = bypass();
      }
    }
    checkModified++;
  }
  destroyRecycleBin(listMemory);
}

bool MetaGraph::trimLongReadTipsSubgraph(){
  bool flagLongRead = false;
  for( IDnum index=1 ; index<=nodeCount(graph) ; index++ ){
    Node* node = getNodeInGraph(graph, index);
    if (node == NULL || subgraphMask[index + nodeCount(graph)] != 1){
      continue;
    }
    if( getUniqueness(node) ){
      continue;
    }
    PassageMarkerI marker, next;
    for( marker=getMarker(node) ; marker!=NULL_IDX ; marker=next ){
      next = getNextInNode(marker);
      flagLongRead = true;
      if( !isInitial(marker) && !isTerminal(marker) ){
	continue;
      }
      if( isTerminal(marker) ){
	marker = getTwinMarker(marker);
      }
      while( !getUniqueness(getNode(marker)) ){
	if( next!=NULL_IDX && (marker==next || marker==getTwinMarker(next)) ){
	  next = getNextInNode(next);
	}
	if( getNextInSequence(marker) != NULL_IDX ){
	  marker = getNextInSequence(marker);
	  destroyPassageMarker( getPreviousInSequence(marker) );
	} else {
	  destroyPassageMarker( marker );
	  break;
	}
      }
    }
  }
  return flagLongRead;
}

void MetaGraph::changeSubgraphMask( int fromID, int toID ){
  for( long i=0 ; i<getSubgraphMaskSize() ; ++i ){
    if( subgraphMask[i] == fromID )
      subgraphMask[i] = toID;
  }
}

long MetaGraph::getSubgraphMaskSize() const {
  return 2 * nodeCount(graph) + 1;
}

bool MetaGraph::checkLongReadExistence() const {
  for( int i=1 ; i<=nodeCount(graph) ; ++i ){
    Node* node = getNodeInGraph(graph, i);
    if( getMarker(node) != NULL_IDX ){
      return true;
    }
  }
  return false;
}

void MetaGraph::resetNodeFlags(){
  resetUniqueness();
  resetNodeStatus( graph );
}

void MetaGraph::resetUniqueness(){
  for( long i=1 ; i<=nodeCount(graph) ; i++ ){
    Node* node = getNodeInGraph( graph, i );
    if( node == NULL )
      continue;
    setUniqueness( node, false );
  }
}

IDnum MetaGraph::getNumActiveNodes() const {
  IDnum count = 0;
  for( IDnum i=1 ; i<=nodeCount(graph) ; i++ ){
    Node* node = getNodeInGraph( graph, i );
    if (node != NULL){
      count++;
    }
  }
  return count;
}
