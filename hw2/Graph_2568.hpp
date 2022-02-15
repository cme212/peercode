#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V, typename E>
class Graph {
 private:
  //vectors that contain unique ID's.
  //When we remove nodes/edges later, we can shift proxy node IDs down
	//but keep these the same.
  //Position in vec will be node position in graph, and value in vector will
	//have unique ID so we can slice into vector of proxynodes/edges
	//(will never delete element from proxynode/edge vecs)
  std::vector<int> nIDconv;//conv for conversion (from proxy id to unique id)
  std::vector<int> eIDconv;
  //Forward declaration of structs used by Node and Edge classes
  struct ProxyNodes;
  struct ProxyEdges;
  struct AdjacentNodes;//used to lessen search in has_edge(a,b))
  //vectors of proxy nodes and proxy edges
  std::vector<ProxyNodes> pnvec;
  std::vector<ProxyEdges> pedgevec;
  //total # of nodes and edges in graph
  //use unsigned to follow example given in class

 public:
  unsigned size_;
  unsigned totalEdges_;

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  typedef V node_value_type;//from stack overflow
  typedef E edge_value_type;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(){
      // total # of nodes and edges are 0 to start off with
      size_=0;
      totalEdges_=0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node:private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
      //already creates invalid node
    }

    /** Return this node's position. */
    const Point& position() const {
      //in this node's graph slice into correct proxy node & return Point att
      //arrow syntax can be used to access elements of a struct:
      //https://www.tutorialspoint.com/What-is-arrow-operator-in-Cplusplus
      return mygraph_->pnvec[nodeIndex_].P;
    }

    //modifiable node position
    Point& position(){
      return mygraph_->pnvec[nodeIndex_].P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //Same idea as position(), but return index attribute
      return mygraph_->pnvec[nodeIndex_].index;
    }

    //return how many nodes are connected to this node.
    size_type degree() const{
      return mygraph_->pnvec[nodeIndex_].ConnectedNodes.size();
    }

    //send nodeIndex of 0
    incident_iterator edge_begin() const{
      return IncidentIterator(nodeIndex_,0,mygraph_);
    }

    //send total # of connected nodes
    incident_iterator edge_end() const{
      return IncidentIterator(nodeIndex_,mygraph_->pnvec[nodeIndex_].\
      ConnectedNodes.size(),mygraph_);
    }

     //return value from ProxyNode struct node
     node_value_type& value(){
      return mygraph_->pnvec[nodeIndex_].value;
    }

    //return value from ProxyNode struct node
    const node_value_type& value() const{
      return mygraph_->pnvec[nodeIndex_].value;
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //if index of this node and node passed in same AND graphs are same
      if(n.nodeIndex_==nodeIndex_ && mygraph_==n.mygraph_){
	return true;
      }else{
         return false;
       }
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      //if this node has index less than index of node passed in
      //will likely decide to change this as we learn more ab how < is used
      if(nodeIndex_<n.nodeIndex_){
        return true;
      }else{
        return false;
       }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //Graph needs a way to construct valid Node objects
    Graph* mygraph_;
    size_type nodeIndex_; //node id for this node
    //since we return Node() instances sometimes (in the add_node method for ex)
    //was throwing errors until I used const cast
    //https://en.cppreference.com/w/cpp/language/const_cast
    Node(const graph_type* graph, size_type index)
	:mygraph_(const_cast<graph_type*>(graph)),nodeIndex_(index) {}

 }; //END OF NODE CLASS

   /** Return the number of nodes in the graph.
    *
    * Complexity: O(1).
    */
    size_type size() const {
      return size_;
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
      return size_;
    }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nodeval = \
  node_value_type()) {
    //initialize vector of adjacent node instances as att of proxynode struct
    std::vector<AdjacentNodes> theseAdj {};
    //add proxynode instance to graph attribute that is vector of proxynodes
    pnvec.push_back(ProxyNodes(position,size_,nodeval,theseAdj));
    //add id to vec that will serve as mapping from proxyid to unique id
    //won't delete from proxynode vec so using size -1 will always be unique
    nIDconv.push_back(pnvec.size()-1);
    (void) position;      // Quiet compiler warning
    size_++;
    //return the newly created node using the private constructor.
    return Node(this, pnvec.size()-1);//use "this" to keep same graph
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //return node with the same id as the one being checked
    Node checkNode = node(n.index());
    //see if these instances are the same
    if(n==checkNode){
      return true;
    }else{
      return false;
     }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if(i<size()){//if index does not exceed size of graph. (calls size method)
      //use Node constructor & vector of ID's that map proxyID to a uniqueID
      return Node(this,nIDconv[i]);
    }else{
      //if the index exceeds the graph size, return an invalid node
      return Node();//invalid node
     }
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      //already invalid edge technically
    }

    edge_value_type& value(){
      return mygraph_->pedgevec[edgeClassIndex_].value;
    }

    const edge_value_type& value() const{
      return mygraph_->pedgevec[edgeClassIndex_].value;
    }

    /** Return a node of this Edge */
    Node node1() const {
      //use node1 attribute of this edge
      return Node(mygraph_,nodeClassIdx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //use node2 attribute of this edge
      return Node(mygraph_,nodeClassIdx2_);
    }

    //return Edge Length
    double length() const{
      Point p1 = node1().position();
      Point p2 = node2().position();
      return norm(p1-p2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //can be equal in either "direction"
      if(node1() == e.node1() && node2() == e.node2()){
         return true;
      }else if(node2() == e.node1() && node1() == e.node2()){
         return true;
      }else{
         return false;
       }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //just like nodes, compare if this edge id less than one we pass in
      if(edgeClassIndex_<e.edgeClassIndex_){
        return true;
      }
      //need comparison if graphs are different
      else if(mygraph_!=e.mygraph_){
        return true;
      }
      //can also compare using the nodes
      else if(node1() < e.node1()){
        return true;
      }
      else{
        return false;
     }
   }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Graph needs a way to construct valid Edge objects
    graph_type* mygraph_;
    //this edge has its own index, and id's for both nodes
    size_type edgeClassIndex_;
    size_type nodeClassIdx1_;
    size_type nodeClassIdx2_;
    //just like in Node constructor, was receiving error until const cast
    Edge(const graph_type* graph,size_type idx,size_type n1,size_type n2)
	:mygraph_(const_cast<graph_type*>(graph)), edgeClassIndex_(idx),\
	 nodeClassIdx1_(n1),nodeClassIdx2_(n2){}

 };//END OF EDGE CLASS

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return totalEdges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    //use vector w/unique id's to define this edge's id
    return Edge(this,eIDconv[i],pedgevec[i].nodeIdx1,pedgevec[i].nodeIdx2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //obtain the index of one node
    size_type  nodeAidx=a.nodeIndex_;
     //each proxy node has vector attribute describing connected nodes
     //only search through this vector for the other node index
     for(unsigned int i = 0; i<pnvec[nodeAidx].ConnectedNodes.size(); i++){
       if(pnvec[nodeAidx].ConnectedNodes[i].connectedNode==b.nodeIndex_)
         return true;
     }
     return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b,const edge_value_type& edgeval = \
  edge_value_type()) {
      //check if edge is already in the graph
      if(has_edge(a,b)==false){
        //add the edge to graph's vector of proxy edges
        pedgevec.push_back(ProxyEdges(totalEdges_,a.nodeIndex_,\
	b.nodeIndex_,edgeval));
        //add this id to the vector of id's that will map to proxy node id's
        eIDconv.push_back(pedgevec.size()-1);
        //for proxy node a, add b as a connected node
        pnvec[a.nodeIndex_].ConnectedNodes.push_back(AdjacentNodes(\
	pedgevec.size()-1,b.nodeIndex_));
        //for proxy node b, add a as a connected node
        pnvec[b.nodeIndex_].ConnectedNodes.push_back(AdjacentNodes(\
	pedgevec.size()-1,a.nodeIndex_));
        totalEdges_++;
        return Edge(this,pedgevec.size()-1,a.nodeIndex_,b.nodeIndex_);
      }
    //if edge exists in graph already, dont inc total num of edges attribute
    //either way, return an Edge instance
    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //reset graph attributes to 0
    size_=0;
    totalEdges_=0;
    //credit to https://www.geeksforgeeks.org/vector-erase-and-clear-in-cpp/
    pedgevec.clear();
    pnvec.clear();
    nIDconv.clear();
    eIDconv.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    //change to forward iterator for std::min_element in shortest path hpp
    using iterator_category = std::forward_iterator_tag;  //Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    //return node with same index as node iterator
    Node operator*() const{
      //use vcgor that maps unique id's to actual
      return Node(itgraph,itgraph->nIDconv[index]);
    }

    //if not already at the end, increase index and return iterator
    NodeIterator& operator++(){
      if(index < itgraph->num_nodes())
      index++;
      return *this;
    }

    //check if all attributes of iterator equal passed in iterator
    bool operator==(const NodeIterator& a) const{
      if(a.index == index && a.itgraph == itgraph)
        return(true);
      else
        return(false);
    }

   //node iterator has its own id and needs to know which graph it refers to
    size_type index;
   private:
    const graph_type* itgraph;
   //allow graph to access attributes
    friend class Graph;
    //constructor
    NodeIterator(size_type i, const graph_type* g): index(i), itgraph(g) {}

  };//END OF NODE ITERATOR CLASS

   //beginning corresponds to index 0
   node_iterator node_begin() const{
     return NodeIterator(0,this);
   }

   //end corresponds to index equal to total size
   node_iterator node_end() const{
     return NodeIterator(this->num_nodes(),this);
   }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    //node 1 from connecting (i.e. incident) edge. node 2 from Adj Node struct
    //id from adjacent node struct
    Edge operator*() const{
      return(Edge(graph,graph->pnvec[beginNode].ConnectedNodes[endNode]\
      .edgeID,beginNode,graph->pnvec[beginNode].ConnectedNodes[endNode]\
      .connectedNode));
    }

    //increase the id and return
    IncidentIterator& operator++(){
      endNode++;
      return(*this);
    }

    //equal if all attributes are equivalent
    bool operator==(const IncidentIterator& a) const{
      if(graph == a.graph && endNode == a.endNode && beginNode == a.beginNode)
        return true;
      else
        return false;
    }

   private:
    friend class Graph;
    //each incident edge connects two nodes
    //also need to know which graph this refers to
    graph_type* graph;
    size_type beginNode;
    size_type endNode;
    //private constructor
    IncidentIterator(size_type begin, size_type end, const graph_type* agraph):
      graph(const_cast <graph_type*>(agraph)), beginNode(begin), endNode(end) {}
 };//END OF INCIDENT ITERATOR CLASS

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    //return edge with id same as iterator id and nodes from proxy edge struct
    Edge operator*() const{
      return(Edge(eitgraph,eitid,eitgraph->pedgevec[eitid].nodeIdx1,\
      eitgraph->pedgevec[eitid].nodeIdx2));
    }

    //increase id and return
    EdgeIterator& operator++(){
      eitid++;
      return *this;
    }

    //equal if all attributes equivalent
    bool operator==(const EdgeIterator& a) const{
      if(eitid==a.eitid && eitgraph==a.eitgraph && proxyEdgeit==proxyEdgeit)
        return(true);
      else
        return(false);
    }

   private:
    friend class Graph;
    //iterator has own id and knows what graph it refers to
    //may also need pointer to proxy edges
    const std::vector<ProxyEdges>* proxyEdgeit;
    size_type eitid;
    graph_type* eitgraph;
    //graph can access attributes
    friend class Graph;
    //private constructor
    EdgeIterator(const graph_type* g, const std::vector<ProxyEdges>* \
    proxyEdgeit2,size_type id2): proxyEdgeit(proxyEdgeit2), eitid(id2),\
    eitgraph(const_cast<graph_type*>(g)) {}

 };//END OF EDGE ITERATOR CLASS

  //beginning if id is 0
  edge_iterator edge_begin() const{
    return (EdgeIterator(this,&pedgevec,0));
  }

  //end if id is same as total edges
  edge_iterator edge_end() const{
    return (EdgeIterator(this,&pedgevec,pedgevec.size()));
  }




  //
  //REMOVING EDGES AND NODES
  //


  /** @brief Remove an edge from the graph and return 1 (if successful)
   * or 0 (if failure).
   * @pre @rnode1 and @rnode2 are valid nodes of current graph
   * @post num_edges() returns old totalEdges_ -1 if successful and returns
   * previous size if failure
   * @post has_edge(rnode1,rnode2) returns false
   * @param [in] rnode1 and rnode2 are valid nodes of graph
   *cannot invalidate edges for pedgevec
   *the old eIDconv entires may not be the same as the proxy Id's
   *Complexity: No more than O(num_nodes()+num_edges())
  */

  size_type remove_edge(const Node& rnode1, const Node& rnode2){

     if(has_edge(rnode1,rnode2)==true){
       //for each of the connected nodes pertaining to node1
       for(unsigned int i=0; i<rnode1.degree();++i){
         //if the connected node is the node in the edge we want to delete
         if(pnvec[rnode1.nodeIndex_].ConnectedNodes[i].connectedNode == \
         rnode2.nodeIndex_){
           //get edge id from proxy struct. this id is unique
           size_type edgeIDproxy=pnvec[rnode1.nodeIndex_].ConnectedNodes[i].\
           edgeID;
           //get actual edge id from proxy edge struct
           size_type edgeUseID=pedgevec[edgeIDproxy].edgeIndex;

           //we need to keep track of the new adjacent nodes
           //delete the adjacent node that relates to the edge we're deleting
           pnvec[rnode1.nodeIndex_].ConnectedNodes.erase(pnvec[rnode1.\
           nodeIndex_].ConnectedNodes.begin()+i);

           // update adjacent nodes from the other node in this edge we delete:
           for (unsigned int i2=0; i2<rnode2.degree(); ++i2){
             if(pnvec[rnode2.nodeIndex_].ConnectedNodes[i2].connectedNode\
              == rnode1.nodeIndex_)
               pnvec[rnode2.nodeIndex_].ConnectedNodes.erase(pnvec[rnode2.\
               nodeIndex_].ConnectedNodes.begin()+i2);
           }

           //update the vector that exists as  map from the unique id to actual
           eIDconv.erase(eIDconv.begin()+edgeUseID);
           //can now still slice into map vec w/ idx & val at position is
             //actual idx

           //shift down idxs of ones that occur after the edge we just deleted
           for(unsigned int i3 = edgeUseID; i3 < eIDconv.size(); ++i3){
             pedgevec[eIDconv[i3]].edgeIndex = i3;
           }

          //update attributes of graph
          totalEdges_ =totalEdges_ -1;
          return 1;
        }//end of "if connected node is one we want to delete"
     }//end of outermost for loop
   }//end of "if we can find the edge in the graph"
    return 0;
  }

  /** @brief Remove an edge from the graph and return 1 (if successful) or
   * 0 (if failure).
   * @pre @redge.node1 and @redge.node2 are valid nodes of current graph
   * @post num_edges() returns old totalEdges_ -1 if successful and returns
   * previous size if failure
   * @post has_edge(rnode1,rnode2) returns false
   * @param[in] = redge is valid edge of graph

   *cannot invalidate edges for pedgevec
   *the old eIDconv entires may not be the same as the proxy Id's
   *Complexity: No more than O(num_nodes()+num_edges())
  */

  //call previously defined function if given an edge
  size_type remove_edge(const Edge& redge){
    //define nodes so we can call the last function
    Node firstN =redge.node1();
    Node secondN =redge.node2();
    size_type successOrFail = remove_edge(firstN, secondN);
    //returns 1 if success and 0 if fail
    return successOrFail;
  }

  /** @brief Remove an edge from the graph and return 1 (if successful) or
   * 0 (if failure).
   * @pre e_it is a valid iterator (dereferencing returns valid edge) of
   * current graph
   * @post num_edges() returns old totalEdges_ -1 if successful and
   * returns previous size if failure
   * @post has_edge(rnode1,rnode2) returns false
   * @param[in] = redge is valid edge iterator of the graph

   *cannot invalidate edges for pedgevec
   *the old eIDconv entires may not be the same as the proxy Id's
   *Complexity: No more than O(num_nodes()+num_edges())
  */

  //call function if given edge iterator
  edge_iterator remove_edge(edge_iterator e_it){
     //define the edge and then call remove_edge with the edge
     //this will then call full remove function
     auto edge = *(e_it);
     remove_edge(edge);
     //re-define edge iterator now that edge was deleted
     edge_iterator continuedIt = e_it;
     //not sure if necessary but makes more clearly how it's updated
     return(continuedIt);
  }

  /** @brief Remove a node from the graph and return 1 (if successful) or
   * 0 (if failure).
   * @pre @rnode is a valid node of current graph
   * @post num_nodes() returns old size_ -1 if successful and returns
   * previous size if failure
   * @post has_node(rnode) returns false
   * @param[in] = rnode is valid node of graph

   *cannot invalidate edges for pnvec
   *the old nIDconv entires may not be the same as the proxy Id's
   *Complexity: No more than O(num_nodes())
  */

  size_type remove_node(const Node& rnode){
    if(has_node(rnode)){
      //per instructions, remove all edges connected to node we want to remove
      while(rnode.degree()>0){
        auto rIncEdge = *(rnode.edge_begin());
        remove_edge(rIncEdge);
      }
      //remove node id from the index conversion vector
      nIDconv.erase(nIDconv.begin() + rnode.index());
      for (size_type i=rnode.index(); i< nIDconv.size(); ++i){
        //shift all the id's down that occured after one that was removed
        pnvec[nIDconv[i]].index=i;
      }

     //just like for edges, update graph info
     size_-=1;
     return 1;
    }
    return 0;
  }

  /** @brief Remove a node from the graph and return 1 (if successful) or
   * 0 (if failure).
   * @pre @rnode is a valid node iterator of current graph
   * @post num_nodes() returns old size_ -1 if successful and returns previous
   * size if failure
   * @post has_node(rnode) returns false
   * @param[in] = n_it is valid node iterator of graph

   *cannot invalidate edges for pnvec
   *the old nIDconv entires may not be the same as the proxy Id's
   *Complexity: No more than O(num_nodes())
  */

  //call previous function if we pass in a node iterator
  node_iterator remove_node(node_iterator n_it){
    //define the node and remove it
    auto node = *n_it;
    remove_node(node);
    //again, not sure if necessary
    //but to make it readable, show how iterator updated
    node_iterator continuedIt = n_it;
    return(continuedIt);
  }



 private:
   //used within an attribute of proxy nodes
   struct AdjacentNodes{
     size_type edgeID;//edge connecting the node
     size_type connectedNode;//adjacent node
     //make constructor for when we add an instance to vec in add_edge
     AdjacentNodes(size_type edgeid, size_type nextNode)
       :edgeID(edgeid),connectedNode(nextNode){}
   };

  //vector of ProxyNodes instances is graph class attribute
  struct ProxyNodes{
    Point P; // struct defined in other .hpp file
    size_type index; //proxy id. will be shifted if nodes removed
    node_value_type value;
    std::vector<AdjacentNodes> ConnectedNodes;//each node has connected nodes
    //make constructor for when we add instance to vector
    //Point can be const as no need to change it
    ProxyNodes(const Point& MyP, size_type myIndex, \
       node_value_type myVal,std::vector<AdjacentNodes> adjNodes)
       :P(MyP),index(myIndex),value(myVal),ConnectedNodes(adjNodes){}
  };

  //use same idea as for nodes since it worked
  struct ProxyEdges{
    //proxy edge id will be shifted if edges removed
    size_type edgeIndex;
    //edges have 2 nodes associated with them
    size_type nodeIdx1;
    size_type nodeIdx2;
    edge_value_type value;
    //make constructor for when we add to vec of proxy edges
    ProxyEdges(size_type eIndex,size_type nIdx1,size_type nIdx2,\
    edge_value_type myVal):edgeIndex(eIndex),nodeIdx1(nIdx1),nodeIdx2(nIdx2),\
    value(myVal){}
  };

};//END OF GRAPH CLASS

#endif // CME212_GRAPH_HPP


