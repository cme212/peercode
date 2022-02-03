#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:
  
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;
  struct adjacent_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Allow the Nodes to support a user-specified value. */
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
  Graph() {
    // HW0: YOUR CODE HERE
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
      graph_ = nullptr;
      nid_ = 0; 
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_-> nodeList[nid_].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value.*/
    node_value_type& value(){
      return graph_ -> nodeList[nid_].nval;
    }

    /** Return a constant reference to this node's value. */
    const node_value_type& value() const{
      return graph_ -> nodeList[nid_].nval;
    }

    /** Return the number of incident edges. */
    size_type degree() const{
      return graph_ -> nodeList[nid_].adj_nodes.size();
    }

    // Start of the incident iterator.
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, nid_, 0);
    }

    // End of incident iterator.
    incident_iterator edge_end () const{
      return IncidentIterator(graph_, nid_, degree()); 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == n.graph_) && (nid_ == n.nid_));
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
      // HW0: YOUR CODE HERE
      assert(graph_ == n.graph_);
      return nid_ < n.nid_;
    }
  

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type nid_;

    Node(const graph_type* graph, size_type nid)
    : graph_(const_cast<graph_type*>(graph)), nid_(nid) {    
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodeList.size();
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
  Node add_node(const Point& position,const node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE
    //Create New node, add its attributes and push to vector
    internal_node new_node;
    new_node.pos = position; 
    new_node.nval = v;
    nodeList.push_back(new_node);

    return Node(this, size() - 1);  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.nid_ < size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);   
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
      graph_ = nullptr; 
      eid_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, firstn_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, secondn_);    
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return ((eid_ == e.eid_) && (graph_ == e.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      assert(graph_ == e.graph_);
      return eid_ < e.eid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_; 
    size_type eid_;
    size_type firstn_;
    size_type secondn_; 

    Edge(const graph_type* graph, size_type eid, size_type firstn, size_type secondn)
    : graph_(const_cast<graph_type*>(graph)), eid_(eid), firstn_(firstn), secondn_(secondn){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edgeList.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type firstn_ = edgeList[i].first;
    size_type secondn_ = edgeList[i].second;
    return Edge(this, i, firstn_, secondn_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //Pull the adjacent node list for node a and check if node b is in the list
    std::vector<adjacent_node> adjNodeList = nodeList[a.nid_].adj_nodes;
    for (unsigned i = 0; i < adjNodeList.size(); ++i){
      if (b.nid_ == adjNodeList[i].adj_nid_){
        return true;
      }
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
    //If edge exists find edge id and return edge
    std::vector<adjacent_node> adjNodeList = nodeList[a.nid_].adj_nodes;
    for (unsigned i = 0; i < adjNodeList.size(); ++i){
      if (b.nid_ == adjNodeList[i].adj_nid_){
        size_type eid_ = adjNodeList[i].iid_;
        return edge(eid_);
      }
    }
    //Push new edge to edge list and create temp adj nodes
    edgeList.push_back(std::make_pair(a.nid_, b.nid_));
    adjacent_node adj1, adj2; 

    //Assign adj node id and edge id and push to the vector of adj nodes
    adj1.adj_nid_ = b.nid_;
    adj1.iid_ = num_edges() - 1;
    adj2.adj_nid_ = a.nid_;
    adj2.iid_ = num_edges() - 1;
    nodeList[a.nid_].adj_nodes.push_back(adj1);
    nodeList[b.nid_].adj_nodes.push_back(adj2);

    return edge(num_edges() - 1);

    }
  

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodeList.clear();
    edgeList.clear();
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
      graph_ = nullptr;
      nid_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the dereferenced node. */
    Node operator*() const {
      return Node(graph_, nid_);
    }

    /** Increment iterator. */
    NodeIterator& operator++(){
      nid_ += 1; 
      return *this;
    }

    /** Test whether this iterator is equal to node iterator. */
    bool operator==(const NodeIterator& n) const{
      return ((graph_ == n.graph_) && (nid_ == n.nid_));
    } 

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type nid_;

    NodeIterator(const graph_type* graph, size_type nid)
    : graph_(const_cast<graph_type*>(graph)), nid_(nid) {    
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return the starting index for node iterator. */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Return the ending index for node iterator. */
  node_iterator node_end() const{
    return NodeIterator(this, num_nodes());
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
      graph_ = nullptr;
      currnid_ = 0;
      adjid_ = 0;

    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return a dereferenced edge. */
    Edge operator*() const{
      //Retrieve edge id and adj node id
      size_type eid_ = graph_ -> nodeList[currnid_].adj_nodes[adjid_].iid_;
      size_type adjn_ = graph_ -> nodeList[currnid_].adj_nodes[adjid_].adj_nid_;
      return Edge(graph_, eid_, currnid_, adjn_);
    }

    /** Increment iterator. */
    IncidentIterator& operator++(){
      adjid_ += 1; 
      return *this;
    }

    /** Test whether this iterator is equal to incident iterator. */
    bool operator==(const IncidentIterator& e) const{
      return ((graph_ == e.graph_) && (currnid_ == e.currnid_) && (adjid_ == e.adjid_));
    } 

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type currnid_;
    size_type adjid_;

    IncidentIterator(const graph_type* graph, size_type currnid, size_type adjid)
    : graph_(const_cast<graph_type*>(graph)), currnid_(currnid), adjid_(adjid) {    
    }
    
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
      graph_ = nullptr;
      eid_ = 0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return the dereferenced edge. */
    Edge operator*() const {
      //Grab node 1 and node 2 id
      size_type nida_ = graph_ -> edgeList[eid_].first;
      size_type nidb_ = graph_ -> edgeList[eid_].second;
      return Edge(graph_, eid_, nida_, nidb_);
    }

    /** Increment iterator. */
    EdgeIterator& operator++(){
      eid_ += 1; 
      return *this;
    }

    /** Test whether this iterator is equal to edge iterator. */
    bool operator==(const EdgeIterator& e) const{
      return ((graph_ == e.graph_) && (eid_ == e.eid_));
    } 


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type eid_;

    EdgeIterator(const graph_type* graph, size_type eid)
    : graph_(const_cast<graph_type*>(graph)), eid_(eid) {    
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return the starting index for edge iterator. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Return the ending index for edge iterator. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  struct internal_node {
    Point pos;
    node_value_type nval;
    std::vector<adjacent_node> adj_nodes; 
  };

  struct adjacent_node{
    size_type adj_nid_;
    size_type iid_;
  };

  //List of internal nodes
  std::vector<internal_node> nodeList;
  //List of edges as a pair of nodes that form the edges
  std::vector<std::pair<size_type, size_type>> edgeList;

};

#endif // CME212_GRAPH_HPP
