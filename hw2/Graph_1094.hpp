#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  typedef V node_value_type;
  typedef E edge_value_type;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {clear();}

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_!=nullptr);
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->nodes_[uid]->position_;
    }
    /** Return modifiable version of node's position*/
    Point& position() {
      assert(graph_!=nullptr);
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->nodes_[uid]->position_;
    }
    /** Return node's initial position **/
    const Point& initialPosition() {
      assert(graph_!=nullptr);
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->nodes_[uid]->initialPosition_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    node_value_type& value(){
      assert(graph_!=nullptr);
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->nodes_[uid]->val_;
    }
    // const node_value_type& value() const;
    const node_value_type& value() const{
      assert(graph_!=nullptr);
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->nodes_[uid]->val_;
    }
    // size_type degree() const;
    size_type degree() const{
      if (graph_==nullptr){
        return 0;
      }
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return graph_->EMap_[uid].size();
    }
    // incident_iterator edge_begin() const;
    IncidentIterator edge_begin() const{
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return IncidentIterator(graph_,nid_,graph_->EMap_[uid].begin());
    }
    // incident_iterator edge_end() const;
    IncidentIterator edge_end() const{
      size_type uid = graph_->nodeMapping_[nid_];
      assert(uid<graph_->nodes_.size());
      return IncidentIterator(graph_,nid_,graph_->EMap_[uid].end());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // if two nodes are invalid, return true
      if (graph_==nullptr && n.graph_==nullptr){
          return true;
      }
      if (graph_==n.graph_ && nid_==n.nid_){
          return true;
      }
      return false;
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
      // if two nodes are equal,return false
      if (*this==n){
          return false;
      }
      // HW0: YOUR CODE HERE
      // if different graph but same idx, return comparison of graph pointers
      if ( nid_ == n.nid_){
        return graph_<n.graph_;
      }
      return nid_<n.nid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to Graph container
    Graph* graph_=nullptr;
    // This node's idx
    size_type nid_;
    // Private Constructor
    Node(const Graph* graph, size_type nid):
        graph_(const_cast<Graph*>(graph)),nid_(nid){}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodeMapping_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
    const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Add new internal_node element to vector storing nodes info
    internal_node* ptr = new internal_node(position,val);
    nodes_.push_back(ptr);
    nodeMapping_.push_back(nodes_.size()-1);
    return Node(this,nodeMapping_.size()-1);        // valid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // Check if in the graph and index within range
    return this==n.graph_ && nodeMapping_.size()>n.nid_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<size());
    return Node(this,i);        // valid node
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *node1_;      // valid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *node2_;      // valid Node
    }

    // Return modifiable value of an edge
    edge_value_type& value(){
      assert(graph_!=nullptr);
      size_type uid = graph_->edgeMapping_[eid_];
      assert(uid<graph_->edges_.size());
      return graph_->edges_[uid]->val_;
    }

    // Return const value of an edge
    const edge_value_type& value() const{
        assert(graph_!=nullptr);
        size_type uid = graph_->edgeMapping_[eid_];
        assert(uid<graph_->edges_.size());
        return graph_->edges_[uid].val_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if in the same graph
      // if two edges are invalid, return true
      if (graph_==nullptr && e.graph_==nullptr){
          return true;
      }
      if (graph_!=e.graph_){
          return false;
      }
      // Check if same nodes in same order
      if (*node1_==*(e.node1_) && *node2_==*(e.node2_)){
          return true;
      }
      // Check if same nodes in reverse order
      if (*node1_==*(e.node2_) && *node2_==*(e.node1_)){
          return true;
      }
      //HW0: YOUR CODE HERE
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //assert(graph_==e.graph_);
      // compare edge index
      if (*this==e){
        return false;
      }
      if (eid_==e.eid_){
        return graph_<e.graph_;
      }
      return eid_<e.eid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph container
    Graph* graph_=nullptr;
    // This is edge's idx
    size_type eid_;
    // Pointer to node 1
    const Node* node1_;
    // Pointer to node 2
    const Node* node2_;
    // Private Constructor
    Edge(const Graph* graph,size_type eid,const Node* node1,const Node* node2):
        graph_(const_cast<Graph*>(graph)),eid_(eid),node1_(node1),node2_(node2){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edgeMapping_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(i<edgeMapping_.size());
    size_type uid = edgeMapping_[i];
    assert(uid<edges_.size());
    return Edge(this,i,edges_[uid]->node1_,edges_[uid]->node2_);        // valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_==this && b.graph_==this);
    auto auid = nodeMapping_[a.nid_];
    auto buid = nodeMapping_[b.nid_];
    if (EMap_.count(auid)>0&&EMap_.at(auid).count(buid)>0){
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_==this && b.graph_==this);
    assert(a.nid_!=b.nid_);
    assert(a.nid_<size()&&b.nid_<size());
    // Check if edge exists
    // Indicator for existance of edge
    bool exist=false;
    // idx of edge if edge existed
    size_type idx = 0;
    // get uid of the nodes
    
    auto auid = nodeMapping_[a.nid_];
    auto buid = nodeMapping_[b.nid_];
    
    // Check if edge exists
    if (EMap_.count(auid)>0&&EMap_.at(auid).count(buid)>0){
        exist = true;
        idx = EMap_.at(auid).at(buid);
    }
    if (exist){
        return Edge(this,idx,&a,&b);;
    }
    // Create new edge if no existing edge
    auto newA = new Node(this,a.nid_);
    auto newB = new Node(this,b.nid_);
    Edge e = Edge(this,edges_.size(),newA,newB);
    // Add length of the edge to edgeVals_
    Point posA = a.position();
    Point posB = b.position();
    auto edge_val = norm(Point(posA.x-posB.x,posA.y-posB.y,posA.z-posB.z));
    edges_.push_back(new internal_edge(newA,newB,edge_val));
    edgeMapping_.push_back(edges_.size()-1);
    // If auid is a key in EMap_
    if (EMap_.count(auid)>0){
        EMap_.at(auid).insert({buid,edgeMapping_.size()-1});
    }
    // If buid is not a key in EMap_
    else {
        EMap_.insert({auid,{{buid,edgeMapping_.size()-1}}});
    }
    // Add (b,a) into EMap_ as well
    // If buid is a key in EMap_
    if (EMap_.count(buid)>0){
        EMap_.at(buid).insert({auid,edgeMapping_.size()-1});
    }
    // If buid_ is not a key in EMap_
    else {
        EMap_.insert({buid,{{auid,edgeMapping_.size()-1}}});
    }

    return e;        // valid Edge
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Delete objects in nodes_
    for (auto i : nodes_){
      if (i!=nullptr){
        delete i;
      }
    }
    nodes_.clear();
    // Delete objects in edges_
    for (auto i : edges_){
      if (i!=nullptr){
        auto n1 = i->node1_;
        auto n2 = i->node2_;
        if (n1!=nullptr){
          delete n1;
        }
        if (n2!=nullptr){
          delete n2;
        }
        delete i;
      }
    }
    edges_.clear();
    EMap_.clear();
    nodeMapping_.clear();
    edgeMapping_.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    Node operator*() const{
      assert(graph_!=nullptr);
      return (*graph_).node(nid_);
    }
    // NodeIterator& operator++()
    NodeIterator& operator++(){
      nid_++;
      return *this;
    }
    // bool operator==(const NodeIterator&) const
    bool operator==(const NodeIterator& nodeIter) const{
      return (graph_==nodeIter.graph_ && nid_==nodeIter.nid_);
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_=nullptr;
    // This node's idx
    size_type nid_;
    // Private Constructor
    NodeIterator(const Graph* graph, size_type nid):
      graph_(const_cast<Graph*>(graph)),nid_(nid){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  NodeIterator node_begin() const{
    return NodeIterator(this,0);
  }
  // node_iterator node_end() const
  NodeIterator node_end() const{
    return NodeIterator(this,num_nodes());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      assert(graph_!=nullptr);
      size_type edgeIdx = iter_->second;
      size_type uid = graph_->edgeMapping_[edgeIdx];
      auto e = *(graph_->edges_[uid]);
      Node n1 = *(e.node1_);
      // if node 1 is our node
      if (n1.index()==nodeId_){
        return (*graph_).add_edge(*(e.node1_),*(e.node2_));
      }
      // if node 1 is not our node, need to flip edge
      else{
        return (*graph_).add_edge(*(e.node2_),*(e.node1_));
      }
    
    }

    // IncidentIterator& operator++()
    IncidentIterator& operator++(){
      iter_++;
      return *this;
    }
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& iter) const{
        return(graph_==iter.graph_&&nodeId_==iter.nodeId_&&iter_==iter.iter_);
    }
   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_=nullptr;
    //idx of the node
    size_type nodeId_;
    //Iterator over EMap[nodeId_]
    std::unordered_map<size_type, size_type>::iterator iter_;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Graph* graph,size_type nid, 
      std::unordered_map<size_type, size_type>::iterator iter):
      graph_(const_cast<Graph*>(graph)),nodeId_(nid),iter_(iter){}
  };

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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      assert(graph_!=nullptr);
      return graph_->edge(eid_);
    }
    // EdgeIterator& operator++()
    EdgeIterator& operator++(){
      eid_++;
      return *this;
    }
    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& iter) const{
      return graph_==iter.graph_ && eid_==iter.eid_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer back to Graph container
    Graph* graph_=nullptr;
    // The edge's idx
    size_type eid_;
    // Private Constructor
    EdgeIterator(const Graph* graph, size_type eid):
      graph_(const_cast<Graph*>(graph)),eid_(eid){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  EdgeIterator edge_begin() const{
    return EdgeIterator(this,0);
  }
  // edge_iterator edge_end() const
  EdgeIterator edge_end() const{
    return EdgeIterator(this,edgeMapping_.size());
  }

  /** Remove an edge from graph, return 1 if edge can be removed, return 0 elsewise
   * @pre @a e is an edge
   * @return 1 if edge can be removed, return 0 otherwise
   * @post has_edge(@a e) == false
   * @post If old has_edge(@a e), new num_edges() == old num_edges()-1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Removed edges will be invalidated
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge& e){
    if(e.eid_>num_edges()){
      return 0;
    }
    auto eid = e.eid_;
    auto n1 = e.node1();
    auto n2 = e.node2();
    // check if is a valid edge
    if (!has_edge(n1,n2)){
      return 0;
    }
    // Remove edge from EMap_
    auto n1uid = nodeMapping_[n1.nid_];
    auto n2uid = nodeMapping_[n2.nid_];
    EMap_[n1uid].erase(n2uid);
    EMap_[n2uid].erase(n1uid);
    
    //swap and pop
    edgeMapping_[eid] = edgeMapping_.back();
    edgeMapping_.pop_back();
  
    // check if edge removed is the last edge in vector
    // update edge idx of moved edge (originally at back) in EMap_
    if (eid<num_edges()){
      auto newE = *(edges_[edgeMapping_[eid]]);
      auto newN1 = *(newE.node1_);
      auto newN2 = *(newE.node2_);
      
      assert(has_node(newN1)&&has_node(newN2));
      auto newN1uid = nodeMapping_[newN1.nid_];
      auto newN2uid = nodeMapping_[newN2.nid_];
  
      EMap_[newN1uid][newN2uid] = eid;
      EMap_[newN2uid][newN1uid] = eid;
    }
    return 1;
  }

  /** Remove an edge from graph given two nodes
   * @pre @a n1 and @a n2 are valid points in the graph
   * @return 1 if edge(n1,n2) can be removed, return 0 otherwise
   * @post has_edge(@a e=(n1,n2)) == false
   * @post If old has_edge(@a e=(n1,n2)), new num_edges() == old num_edges()-1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Removed edges will be invalidated, nodes still remain valid
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& n1, const Node& n2){
    if (!has_edge(n1,n2)){
      return 0;
    }
    auto e = add_edge(n1,n2);
    assert(e.eid_<num_edges());
    return remove_edge(e);
  }


  /** Remove an edge from graph given an edge_iterator
   * @pre @a e_it is an edge_iterator
   * @return a valid edge iterator if *e_it can be removed
   * return invalid interator otherwise
   * @post has_edge(@a *e_it) == false
   * @post If old has_edge(@a *e_it), new num_edges() == old num_edges()-1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Edge iterator will be invalidated, corresponding edge will be invalidated
   *
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it){
    if (e_it==edge_end()){
      return e_it;
    }
    auto result = remove_edge(*e_it);
    if (result==0){
      return edge_end();
    }
    return edge_begin();
  }

  /** Remove a node from the graph
   * @pre @a n is a node
   * @return 1 if node can be removed, 0 otherwise
   * @post has_node(@a n) == false, has_edge(edge incident to n)==false
   * @post If old has_node(@a n), new size() == old size()-1
   *       Else,                        new size() == old size()
   *
   * Invalidate node n  and edges incident to the node
   *
   * Complexity: O(num nodes())
   */
  size_type remove_node(const Node& n){
    if(!has_node(n)){
      return 0;
    }
    auto nid = n.index();
    for (auto iter = n.edge_begin();iter!=n.edge_end();){
      auto e = *iter;
      ++iter;
      remove_edge(e);
    }
    //swap
    nodeMapping_[nid] = nodeMapping_.back();
    
    // check if node removed is the last node in nodeMapping
    // change node value of moved node(originally at back) in internal_edge
    if (nid<size()-1){
      for (auto iter = n.edge_begin();iter!=n.edge_end();++iter){
        auto e = *iter;
        auto eid = e.eid_;
        auto uid = edgeMapping_[eid];
        auto int_edge = edges_[uid];
        if (int_edge->node1_->index()==nodeMapping_.size()-1){
          delete int_edge->node1_;
          int_edge->node1_ = new Node(this,nid);
        }
        if (int_edge->node2_->index()==nodeMapping_.size()-1){
          delete int_edge->node2_;
          int_edge->node2_ = new Node(this,nid);
        }
      }
      
    }
    // pop
    nodeMapping_.pop_back();
    return 1;
  }

  /** Remove a node from the graph
   * @pre @a n_it is a node iterator
   * @return valid iterator if node can be removed, node_end() otherwise
   * @post has_node(@a *n_it) == false,has_edge(edge incident to n)==false
   * @post If old has_node(@a *n_it), new size() == old size()-1
   *       Else,                        new size() == old size()
   *
   * Invalidate the node iterator, corresponding node and edges incident to the node
   *
   * Complexity: O(num nodes())
   */
  node_iterator remove_node(node_iterator n_it){
    if (n_it==node_end()){
      return n_it;
    }
    auto result = remove_node(*n_it);
    if (result==0){
      return node_end();
    }
    return node_begin();
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //Internal type for node elements
  struct internal_node {
      Point position_;
      Point initialPosition_;
      node_value_type val_;
      internal_node(const Point& position, const node_value_type& val){
          position_ = Point(position.x,position.y,position.z);
          initialPosition_ = Point(position.x,position.y,position.z);
          val_=val;
      }
  };

  struct internal_edge {
    const Node* node1_;
    const Node* node2_;
    edge_value_type val_;

    internal_edge(const Node* node1, const Node* node2, const edge_value_type& val){
      node1_ = node1;
      node2_ = node2;
      val_ = val;
    }
  };
  // vector for storing info of nodes
  std::vector<internal_node*> nodes_;

  // For quick access of edge by index
  std::vector<internal_edge*> edges_;
  // For quick check for existence of edge
  // keys are uid of nodes and value is idx of edge
  std::unordered_map<size_type,std::unordered_map<size_type,size_type>> EMap_;

  //Stores mapping of active nodes
  std::vector<size_type> nodeMapping_;
  //Stores mapping of active edges
  std::vector<size_type> edgeMapping_;

};

#endif // CME212_GRAPH_HPP