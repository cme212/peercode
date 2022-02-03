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
template <typename V = int>
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
  using node_value_type = V;
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
  Graph(): nsize_(0), next_nid_(0),esize_(0),next_eid_(0) {
    //std::cout<<"Size of Node: "<<sizeof(Node)<<std::endl;
    //std::cout<<"Size of Edge: "<<sizeof(Edge)<<std::endl;
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
      return graph_->nodes_[nid_]->position_;
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
      return graph_->nodes_[nid_]->val_;
    }
    // const node_value_type& value() const;
    const node_value_type& value() const{
      assert(graph_!=nullptr);
      return graph_->nodes_[nid_]->val_;
    }
    // size_type degree() const;
    size_type degree() const{
      if (graph_==nullptr){
        return 0;
      }
      return graph_->nodes_[nid_]->incidentEIndex_.size();
    }
    // incident_iterator edge_begin() const;
    IncidentIterator edge_begin() const{
        return IncidentIterator(graph_,nid_,0);
    }
    // incident_iterator edge_end() const;
    IncidentIterator edge_end() const{
        return IncidentIterator(graph_,nid_,degree());
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
    // This node's unique identification number
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
    return nsize_;
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
    ++nsize_;
    ++next_nid_;
    return Node(this,next_nid_-1);        // valid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // Check if in the graph and index within range
    return this==n.graph_ && nsize_>n.nid_;
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
      assert(graph_==e.graph_);
      // compare edge index
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
    // This is edge's unique identification number
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
    return esize_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(i<esize_);
    return *edges_[i];        // valid Edge
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
    if (EMap_.count(a.nid_)>0&&EMap_.at(a.nid_).count(b.nid_)>0){
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
    // Check if edge exists
    // Indicator for existance of edge
    bool exist=false;
    // idx of edge if edge existed
    size_type idx = 0;
    // Check if edge exists
    if (EMap_.count(a.nid_)>0&&EMap_.at(a.nid_).count(b.nid_)>0){
        exist = true;
        idx = EMap_.at(a.nid_).at(b.nid_);
    }
    if (exist){
        return Edge(this,idx,&a,&b);;
    }
    // Create new edge if no existing edge
    Edge* e = new Edge(this,next_eid_,&a,&b);
    edges_.push_back(e);
    // If a.nid_ is a key in EMap_
    if (EMap_.count(a.nid_)>0){
        EMap_.at(a.nid_).insert({b.nid_,next_eid_});
    }
    // If a.nid_ is not a key in EMap_
    else {
        EMap_.insert({a.nid_,{{b.nid_,next_eid_}}});
    }
    // Add (b,a) into EMap_ as well
    // If b.nid_ is a key in EMap_
    if (EMap_.count(b.nid_)>0){
        EMap_.at(b.nid_).insert({a.nid_,next_eid_});
    }
    // If b.nid_ is not a key in EMap_
    else {
        EMap_.insert({b.nid_,{{a.nid_,next_eid_}}});
    }
    // Add edge index to incident edge vector of internal_node
    nodes_[a.nid_]->incidentEIndex_.push_back(next_eid_);
    nodes_[b.nid_]->incidentEIndex_.push_back(next_eid_);
    ++esize_;
    ++next_eid_;
    return *e;        // valid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Delete objects in nodes_
    for (auto& i : nodes_){
        delete i;
    }
    nodes_.clear();
    // Delete objects in edges_
    for (auto& i : edges_){
        delete i;
    }
    edges_.clear();
    EMap_.clear();
    nsize_=0;
    next_nid_=0;
    esize_=0;
    next_eid_=0;
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
    // This node's unique identification number
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
      size_type edgeIdx = graph_->nodes_[nodeId_]->incidentEIndex_[index_];
      Edge e = *(graph_->edges_[edgeIdx]);
      Node n1 = e.node1();
      // if node 1 is our node
      if (n1.index()==nodeId_){
        return e;
      }
      // if node 1 is not our node, need to flip edge
      else{
        return (*graph_).add_edge(e.node2(),e.node1());
      }
    }
    // IncidentIterator& operator++()
    IncidentIterator& operator++(){
      index_++;
      return *this;
    }
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& iter) const{
        return(graph_==iter.graph_&&nodeId_==iter.nodeId_&&index_==iter.index_);
    }
   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_=nullptr;
    //Index of the node
    size_type nodeId_;
    //index of position in incidentEIndex_ of the internal_node
    size_type index_;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Graph* graph,size_type nid, size_type index):
      graph_(const_cast<Graph*>(graph)),nodeId_(nid),index_(index){}
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
      return *(graph_->edges_[eid_]);
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
    // The edge's unique identification number
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
    return EdgeIterator(this,esize_);
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //Internal type for node elements
  struct internal_node {
      // position of node
      Point position_;
      // value of node
      node_value_type val_;
      //index of edges incident to node
      std::vector<size_type> incidentEIndex_;
      //constructor
      internal_node(const Point& position, const node_value_type& val){
          position_ = Point(position.x,position.y,position.z);
          val_=val;
      }
  };
  // vector for storing info of nodes
  std::vector<internal_node*> nodes_;
  // number of nodes
  size_type nsize_;
  size_type next_nid_;
  // For quick access of edge by index
  std::vector<Edge*> edges_;
  // For quick check for existence of edge
  std::unordered_map<size_type,std::unordered_map<size_type,size_type>> EMap_;
  // number of edges
  size_type esize_;
  size_type next_eid_;

};

#endif // CME212_GRAPH_HPP
