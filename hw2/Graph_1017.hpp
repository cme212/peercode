#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
 private:
  //declare internal types for Graph
  struct internal_node;
  struct internal_edge;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  typedef V node_value_type;
  
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  typedef E edge_value_type;
  
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
  using map_type = std::unordered_map<size_type, size_type>;
  using set_type = std::set<size_type>;
  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), incident_map_(),
      active_nodes_(), active_edges_(),
      next_node_id_(0), next_edge_id_(0) {}

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
    Node() : graph_(nullptr) {}

    /** Return this node's position. */
    const Point& position() const {
      //if invalid node, no position
      assert(graph_);
      //find the corresponding internal node in nodes_
      return graph_->nodes_[node_id_].point;
    }
    
    /** Return this node's position */
    Point& position() {
      //if invalid node, no position
      assert(graph_);
      //find the corresponding internal node in nodes_
      return graph_->nodes_[node_id_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //if invalid node, no index
      assert(graph_);
      //find iterator pointing to node_id_ in active_nodes_
      auto it = graph_->active_nodes_.find(node_id_);
      size_type idx = 0;
      while (it != graph_->active_nodes_.begin()) {
        it--; idx++;
      }
      assert(idx < graph_->size());
      return idx; //return index
    }

    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the value of the node. if invalid, raises assert error.
     */
    node_value_type& value() {
      //if invalid node, no value
      assert(graph_);
      assert(node_id_ < graph_->nodes_.size());
      return graph_->nodes_[node_id_].nval;
    }
    
    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the value of the node as a const.
     *         If invalid, raises assert error.
     */
    const node_value_type& value() const {
      //if invalid node, no value
      assert(graph_);
      assert(node_id_ < graph_->nodes_.size());
      return graph_->nodes_[node_id_].nval;
    };
    
    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the number of edges for this node size.
     * @post 0 <= return <= num_nodes() - 1
     */
    size_type degree() const {
      //if invalid node or valid node doesn't have any edges
      if (!graph_ or
          (graph_->incident_map_.find(node_id_)==graph_->incident_map_.end())){
        return 0;
      }
      //if node is in incident_map
      const auto& m = graph_->incident_map_.at(node_id_);
      return m.size();
    }
    
    /** Return an incident iterator pointing to the first edge of the node.
     * @pre current node is a valid node in this graph.
     * @return an incident iterator pointing to this node's
     * first edge. Return is an invalid incident iterator,
     * if a node doesn't have any edges.
     */
    IncidentIterator edge_begin() const {
      //if invalid node, then no edges nor incident iterator
      assert(graph_);
      //this node isn't in incident_map_; return an invalid IncidentIterator
      //     but inc.iterator's graph ptr would be nullptr?
      if (graph_->incident_map_.find(node_id_)==graph_->incident_map_.end()) {
        return IncidentIterator();
      }
      //if node is in incident_map
      const auto& m = graph_->incident_map_.at(node_id_);
      return IncidentIterator(graph_, node_id_, m.begin());
    }
    /** Return an incident iterator pointing past the last edge of the node.
     * @pre current node is a valid node in this graph.
     * @return an incident iterator pointing to this graph's this node
     *         and past its last edge.
     * @return return is an invalid incident iterator if a node doesn't have
     *         any edges.
     */
    IncidentIterator edge_end() const {
      //if invalid node, then no edges nor incident iterator
      assert(graph_);
      //this node isn't in incident_map_; return an invalid IncidentIterator
      if (graph_->incident_map_.find(node_id_)==graph_->incident_map_.end()) {
        return IncidentIterator();
      }
      //if node is in incident_map
      const auto& m = graph_->incident_map_.at(node_id_);
      return IncidentIterator(graph_, node_id_, m.end());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (!graph_ and !n.graph_) {
        return true;
      }
      // return true if same graph, same index
      return (graph_ == n.graph_ and node_id_ == n.node_id_);
    }
    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy:
     * For any two nodes x and y, exactly one of x == y,
     * x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (!graph_ and !n.graph_) {
        return false;
      }
      // Comparison based on node_id_
      if (node_id_ != n.node_id_) {
        return node_id_ < n.node_id_;
      }
      return graph_ < n.graph_;
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //Pointer to the Graph this Node belongs to
    Graph* graph_;
    //This Node's unique identifier
    size_type node_id_;
    //Private Constructor
    Node(const Graph* graph, size_type node_id)
         : graph_(const_cast<Graph*> (graph)), node_id_(node_id) {}
  };
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // size of (active) nodes_ in this graph
    return active_nodes_.size();
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
  Node add_node (const Point& position,
                 const node_value_type& val = node_value_type()) {

    internal_node intNode {position, next_node_id_, val};
    nodes_.push_back(intNode);
    active_nodes_.insert(next_node_id_);
    //add to incident map with val=empty map
    if (incident_map_.find(next_node_id_) == incident_map_.end()) {
      map_type empty_map;
      incident_map_.insert(std::make_pair(next_node_id_, empty_map));
    }
    
    ++next_node_id_;
    return Node(this, next_node_id_-1);
  }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_; //check if n points to this Graph
  }
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    //i is index in active_nodes_ -> get nodeID
    auto node_it = active_nodes_.begin();
    while (i > 0) {
      node_it++; i--;
    }
    return Node(this, *node_it);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes.
   * Two Edges with the same nodes are considered equal,
   * if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr) {}

    /** Return a node of this Edge */
    Node node1() const {
      assert(graph_);
      return Node(graph_, first_node_id_);
    }
    /** Return the other node of this Edge */
    Node node2() const {
      assert(graph_);
      return Node(graph_, second_node_id_);
    }
    /** Return the l2-norm of this edge's length*/
    double length() const {
      if (!graph_) {
        return 0;
      }
      return norm_2(node1().position() - node2().position());
    }
    /**Return the edge_value_type value of this edge, possibly the length*/
    edge_value_type& value() {
      assert(graph_);
      assert(edge_id_ < graph_->edges_.size());
      return graph_->edges_[edge_id_].eval;
    }
    /**Return const edge_value_type value of this edge, possibly the length*/
    const edge_value_type& value() const {
      assert(graph_);
      assert(edge_id_ < graph_->edges_.size());
      return graph_->edges_[edge_id_].eval;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (!graph_ and !e.graph_) {
        return true;
      }
      // same two nodes
      return graph_ == e.graph_ and
              ((node1()==e.node1() && node2()==e.node2()) or
              (node1()==e.node2() && node2()==e.node1()));
    }
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (!graph_ and !e.graph_) {
        return false;
      }
      // Compared based on edge_id_
      if (edge_id_ != e.edge_id_) {
        return edge_id_ < e.edge_id_;
      }
      return graph_ < e.graph_;
    }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //This edge's unique identifier
    size_type edge_id_;
    //This edge's nodes specified in order.
    size_type first_node_id_;
    size_type second_node_id_;
    //Private Constructor
    Edge(const Graph* graph,
         size_type edge_id,
         size_type first_node_id,
         size_type second_node_id)
         : graph_(const_cast<Graph*> (graph)),
           edge_id_(edge_id),
           first_node_id_(first_node_id),
           second_node_id_(second_node_id) {}
  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Number of active edges in this graph
    return active_edges_.size();
  }
  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    //get the index for this Edge @ idx=i
    auto edge_it = active_edges_.begin();
    while (i > 0) {
      edge_it++; i--;
    }
    size_type eid = *edge_it;
    internal_edge intEdge = edges_[eid];
    return Edge(this, eid, intEdge.nodeA_id, intEdge.nodeB_id);
  }
  
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.graph_ and b.graph_); //check a and b are valid Nodes
    //if a is in incident_map_, b must be in it, too.
    const auto it = incident_map_.find(a.node_id_);
    if (it != incident_map_.end()) {
      const auto& map_a = it->second;
      return map_a.find(b.node_id_) != map_a.end();
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else, new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    //CASE A: Edge{a,b} already exists
    if (has_edge(a, b)) {
      const auto& map_a = incident_map_.at(a.node_id_);
      const auto it = map_a.find(b.node_id_);
      size_type eid = it->second;
      return Edge(this, eid, a.node_id_, b.node_id_);
    }
    //CASE B: New Edge
    //add edge to edges_ vector
    internal_edge intEdge {a.node_id_, b.node_id_, next_edge_id_,
                           static_cast<edge_value_type>
                           (norm_2(a.position() - b.position())) };
    edges_.push_back(intEdge);
    //add edgeID to active_edges_
    active_edges_.insert(next_edge_id_);
    
    //Assumes Nodes A and B are already in this graph (!!)
    //1. insert Node B to Node A's incident map
    auto& a_map = incident_map_.at(a.node_id_);
    a_map.insert(std::make_pair(b.node_id_, next_edge_id_));
    //2. insert Node A to Node B's incident map
    auto& b_map = incident_map_.at(b.node_id_);
    b_map.insert(std::make_pair(a.node_id_, next_edge_id_));
    
    ++next_edge_id_;
    return Edge(this, next_edge_id_-1, a.node_id_, b.node_id_);
  }
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset all private members of this Graph
    next_node_id_ = 0;
    next_edge_id_ = 0;
    nodes_.clear();
    edges_.clear();
    active_nodes_.clear();
    active_edges_.clear();
    incident_map_.clear();
  }

  //
  // Node Iterator
  //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered <NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  //Weak Category, Proxy
    typedef typename set_type::const_iterator set_iterator_;
    
    /** Construct an invalid NodeIterator. */
    NodeIterator(): graph_(nullptr) {}

    /** Return a node pointed to by this NodeIterator.
     * @pre node iterator points to a valid internal node in nodes_ .
     * @return a Node with the node_id specified by the internal node.
     */
    Node operator*() const {
      //if this iterator is invalid, raise error.
      assert(graph_);
      return Node(graph_, *nit_);
    }
    /** Return this NodeIterator by reference after incrementing.
     * @pre node iterator points to a valid internal node in nodes_
     *      in this graph, not past it.
     * @return this NodeIterator by reference, incremented from pre-call.
     */
    NodeIterator& operator++() {
      //if this iterator is invalid, raise error.
      assert(graph_);
      nit_++;
      return *this;
    }
    /** Return bool indicating if this NodeIterator equals other.
     * @pre NodeIterators are valid.
     * @return 1 if they point to the same graph and element in nodes_,
     *         or  if they're both invalid iterators.
     *         0 otherwise.
     */
    bool operator==(const NodeIterator& other) const {
      if (!graph_ and !other.graph_) {
        return true;
      }
      return nit_ == other.nit_ and graph_ == other.graph_;
    }
   private:
    friend class Graph;
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //Iterator over nodes_ vector for graph_
    set_iterator_ nit_;
    NodeIterator(const Graph* graph, const set_iterator_ nit)
                 : graph_(const_cast<Graph*> (graph)), nit_(nit) {}
  };

  /** Return a NodeIterator pointing to the first elem of nodes_.
   * @pre nodes_ is a valid vector. its size must be >= 0.
   * @return a NodeIterator pointing to the beginning of nodes_.
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, active_nodes_.begin());
  }
  /** Return a NodeIterator pointing past the last elem of nodes_.
   * @pre nodes_ is a valid vector. its size must be >= 0.
   * @return a NodeIterator pointing to the end of nodes_.
   */
  NodeIterator node_end() const {
    return NodeIterator(this, active_nodes_.end());
  }

  //
  // Incident Iterator
  //
  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    typedef typename map_type::const_iterator map_iterator_ ;
    
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr) {}

    /** Return an Edge pointed to by this IncidentIterator.
     * @pre edge iterator points to a valid internal edge in this graph.
     * @pre egde iterator is valid.
     * @return an Edge with the edge ID from the internal edge it pts to.
     * @post returned edge must return this nodeID_ as its 1st node.
     */
    Edge operator*() const {
      assert(graph_); //confirm this iterator is valid.
      
      if (pairs_ == graph_->incident_map_.at(nodeID_).end()) {
        return Edge();
      }
      size_type second_nodeID_ = (*pairs_).first;
      size_type eid = (*pairs_).second;
      return Edge(graph_, eid, nodeID_, second_nodeID_);
    }
    /** Return this after incrementing.
     * @pre edge iterator is valid.
     * @return return == this with pairs_ incremented.
     */
    IncidentIterator& operator++() {
      assert(graph_); //confirm this iterator is valid.
      pairs_++;
      return *this;
    }
    /** Return a boolean indicating if this and other iterators are equal.
     * @return 1 if they point to the same graph and element in edges_,
     *         or if they're both invalid iterators.
     *         0 otherwise.
     */
    bool operator==(const IncidentIterator& other) const {
      if (!graph_ and !other.graph_) {
        return true;
      }
      return graph_ == other.graph_
             && nodeID_ == other.nodeID_
             && pairs_ == other.pairs_;
    }
   private:
    friend class Graph;
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //nodeID whose edges this iterator goes over.
    int nodeID_;
    //iterator for the edges.
    map_iterator_ pairs_;
    
    IncidentIterator(const Graph* graph, int nID, const map_iterator_ p)
            : graph_(const_cast<Graph*> (graph)), nodeID_(nID), pairs_(p) {}
  };

  //
  // Edge Iterator
  //
  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered <EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    typedef typename set_type::const_iterator set_iterator_;
    
    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr) {}

    /** Return an edge this is pointing to.
     * @pre this is a valid iterator.
     * @return an edge pointed to by this iterator.
     * @post return==a valid edge.
     */
    Edge operator*() const {
      assert(graph_); //confirm that this is a valid iterator.
      
      internal_edge intEdge = graph_->edges_[*eit_];
      return Edge(graph_, intEdge.edge_id,
                  intEdge.nodeA_id, intEdge.nodeB_id);
    }
    /** Return this.
     * @pre this is a valid iterator.
     * @return this after incrementing @a eit_.
     * @post return is a valid iterator.
     */
    EdgeIterator& operator++() {
      assert(graph_); //confirm that this is a valid iterator.
      eit_++;
      return *this;
    }
    /** Return a boolean indicating if 2 iterators are equal..
     * @return 1 if 2 iterators' graph_ and eit_  are equal.
     * @return 0 if they're not.
     * @return if both this and other are invalid, return 1.
     */
    bool operator==(const EdgeIterator& other) const {
      //if both iterators are invalid, return true
      if (!graph_ and !other.graph_) {
        return true;
      }
      return graph_ == other.graph_ and eit_ == other.eit_;
    }

   private:
    friend class Graph;
    //Pointer to the Graph this Node belongs to
    Graph* graph_;
    //Iterator pointing to an int. edge in edges_
    set_iterator_ eit_;
    
    EdgeIterator(const Graph* graph, const set_iterator_ eit)
                 : graph_(const_cast<Graph*> (graph)), eit_(eit) {}
  };
  /** Return an EdgeIterator pointing at the 1st elem of edges_.
   * @pre active_edges_ is a valid vector. its size must be >= 0.
   * @return an EdgeIterator pointing to the 1st elem of active_edges_.
   * @post a valid EdgeIterator
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, active_edges_.begin());
  }
  /** Return an EdgeIterator pointing at the last elem of edges_.
   * @pre active_edges_ is a valid vector. its size must be >= 0.
   * @return an EdgeIterator pointing to the last elem of active_edges_.
   * @post a valid EdgeIterator
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, active_edges_.end());
  }
  
  /** Return a boolean value indicating if the removal of @a n was succesful.
   * @param[in] n: node being deleted. @a n's value itself is not modified.
   * @pre Node @a n is a valid Node in this Graph.
   * @return Returns 1 if @a n is successfully removed from this graph.
   *         Else, returns 0, when n is not an active Node.
   * @post @a n is removed from this graph and active_nodes_.
   *       All edges connected to @a n are removed.
   * @post num_nodes() = old num_nodes() - 1, if successful.
   * Complexity: O(n's degree())
   *          Erasing @a n from active_nodes_ takes linear time,
   *          which is longer than iterating over @a n's edges.
   */
  size_type remove_node(const Node& n) {
    if (n.graph_) { //check n is valid
      auto nit = active_nodes_.find(n.node_id_);
      //if n is not an active Node, return 0
      if (nit == active_nodes_.end()) {
        return 0;
      }
      //remove all edges to n
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        size_type n2 = (*it).node2().node_id_; //nodeID of the other node
        //Erase n from n2's incident map
        incident_map_[n2].erase(n.node_id_);
        //Remove edge from active edges
        active_edges_.erase((*it).edge_id_);
      }
      incident_map_.erase(n.node_id_);
      active_nodes_.erase(n.node_id_);
      return 1;
    }
    return 0;
  }
  
  /** Return a NodeIterator pointing to the 1st active Node of this graph.
   * @param[in] n_it: valid NodeIterator pointing to a node being deleted.
   * @pre Node @a n_it is a valid NodeIterator.
   * @return Returns @a n_it pointing to the 1st elem in active_nodes_,
   *         if removal is successful. Else, returns an invalid NodeIterator.
   * @post @a *n_it  is removed from this graph and active_nodes_.
   *       All edges connected @a *n_it it are removed.
   * @post num_nodes() = old num_nodes() - 1, if successful.
   * Complexity: O(degree() )
   *          calls size_type remove_node(const Node&).
   */
  NodeIterator remove_node(NodeIterator n_it) {
    assert(n_it.graph_); //check NodeIterator is valid
    if (n_it != node_end()) {
      remove_node(*n_it);
    }
    return node_begin();
  }
  /** Return a boolean indicating if the removal of Edge connecting
   *  @a a,b  was succesful.
   * @param[in] a,b: valid Nodes whose edge  is being deleted.
   * @pre Nodes @a a,b are valid Nodes and have a connecting edge.
   * @return Returns 1 if @a Edge(a,b) is successfully removed from this graph.
   *         Else, returns 0 - if Edge is not active in the graph.
   * @post @a Edge(a,b)  is removed from this graph, active_edges_, and
   *       incident_map_. Node @a a is removed from Node @a b 's
   *       incident map, and vice versa.
   * @post num_edges() = old num_edges() - 1, if successful.
   * Complexity: O(1); erasing a single element from a vector
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (a.graph_ and b.graph_ and has_edge(a,b)) {
      //Get edgeID before removing edges
      size_type eid = incident_map_[a.node_id_][b.node_id_];
      //Erase a and b from its respective incident map
      int erase1 = incident_map_.at(a.node_id_).erase(b.node_id_);
      int erase2 = incident_map_.at(b.node_id_).erase(a.node_id_);
      int erase3 = active_edges_.erase(eid);//Remove edge from active edges
      //erase not successful
      if (erase1==0 or erase2==0 or erase3==0) {
        return 0;
      }
      
      return 1;
    }
    return 0;
  }
  /** Return a boolean indicating if the removal of Edge @a e was succesful.
   * @param[in] e: Edge being deleted.
   * @pre Edge @a e is a valid and active Edge in this graph.
   * @return Returns 1 if @a e is successfully removed from this graph.
   *         Else, returns 0.
   * @post @a e  is removed from this graph, active_edges_, and
   *       incident_map_, by calling remove_edge on its 2 Nodes.
   * Complexity: O(1)
   *          calls size_type remove_edge(const Node&, const Node&)
   */
  size_type remove_edge(const Edge& e) {
    assert(e.graph_); //check e is a valid Edge
    return remove_edge(e.node1(), e.node2());
  }
  /** Return an EdgeIterator pointing to the 1st active Edge of this graph.
   * @param[in] e_it: EdgeIterator pointing to an edge being deleted.
   * @pre @a e_it is a valid EdgeIterator pointing at an active Edge.
   * @return Returns @a e_it pointing at the 1st edge in active_edges_,
   *         if edge is successfully removed.
   *         Else, returns an invalid EdgeIterator.
   * @post @a *e_it  is removed from this graph, active_edges_,
   *       and incident_map_.
   * @post num_edges() = old num_edges() - 1, if successful.
   * Complexity: O(1)
   *          calls size_type remove_edge(const Node&, const Node&)
   */
  EdgeIterator remove_edge(EdgeIterator e_it) {
    assert(e_it.graph_); //check EdgeIterator is valid
    if (e_it != edge_end()) {
      remove_edge(*e_it);
    }
    return edge_begin();
  }
  
private:
  struct internal_node {
    Point point; //the Point the node stores
    size_type node_id; //unique identification of the node
    node_value_type nval;
  };
  struct internal_edge {
    size_type nodeA_id;
    size_type nodeB_id;
    size_type edge_id; //unique identification of the edge
    edge_value_type eval;
  };
  
  //stores all nodes as internal nodes
  std::vector<internal_node> nodes_;
  //stores all edges as internal edges
  std::vector<internal_edge> edges_;
  //stores key=node, val=map of (node2, edgeID) pairs
  std::unordered_map<size_type, map_type> incident_map_;
  //stores node IDs of active nodes in this graph
  set_type active_nodes_;
  //stores edge IDs of active edges in this graph
  set_type active_edges_;
  
  size_type next_node_id_;
  size_type next_edge_id_;
};

#endif // CME212_GRAPH_HPP
