#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Synonym for V. */
  using node_value_type = V;

  /** Synonym for E. */
  using edge_value_type = E;

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
     * is occasionally useful to declare an @i invalid node, and assign
     * a valid node to it later. For example:
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
    Node() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point &position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].position_;
    }

    /** Return this node's position; non-const. */
    Point &position() { return graph_->nodes_[uid_].position_; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value, of type V. */
    node_value_type &value() { return graph_->nodes_[uid_].value_; }

    /** Return this node's value, of type V; const version */
    const node_value_type &value() const { return graph_->nodes_[uid_].value_; }

    /** Return number of edges connected to this node */
    size_type degree() const {
      // Check if a Node has edges incident to it
      bool connected = graph_->mapping_.find(uid_) != graph_->mapping_.end();
      return connected ? graph_->mapping_.at(uid_).size() : 0;
    }

    /** Returns IncidentIterator at pointing to the first element of
     * incident edges to this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->mapping_[uid_].begin(),
                              graph_->mapping_[uid_].end());
    }
    /** Returns IncidentIterator referring to the past-the-end element
     * of incident edges to this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->mapping_[uid_].end(),
                              graph_->mapping_[uid_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two
     * nodes x and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_) {
        return uid_ < n.uid_;
      }
      return graph_ < n.graph_;
    }

    /** Test whether this node is valid. */
    bool valid() const {
      return uid_ >= 0 && uid_ < nodes_.size() &&
             nodes_[uid_].idx_ < i2u_nodes_.size() &&
             i2u_nodes_[nodes_[uid_].idx_] == uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for
    // Node that will not be visible to users, but may be useful within
    // Graph. i.e. Graph needs a way to construct valid Node objects
    Graph *graph_;   // Pointer back to the Graph
    size_type uid_;  // This Node's index in nodes_

    /** Private Constructor for a valid node. */
    Node(const Graph *graph, const size_type &uid)
        : graph_(const_cast<Graph *>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return i2u_nodes_.size(); }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type &value = node_value_type()) {
    size_type uid = nodes_.size();      // New node's uid
    size_type idx = i2u_nodes_.size();  // New node's idx
    i2u_nodes_.push_back(uid);
    nodes_.push_back(internal_node(position, value, idx));
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_ && is_valid_node(n.uid_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, i2u_nodes_[i]); }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same
   * nodes are considered equal if they connect the same nodes, in either
   * order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_uid_);
    }

    /** Return the Euclidean distance between the Edge's nodes */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return this edge's value, of type E. */
    edge_value_type &value() { return graph_->edges_[uid_].value_; }

    /** Return this edges's value, of type E; const version */
    const edge_value_type &value() const { return graph_->edges_[uid_].value_; }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      // HW0: YOUR CODE HERE
      return graph_ == e.graph_ && uid_ == e.uid_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // HW0: YOUR CODE HERE
      if (graph_ == e.graph_) {
        return uid_ < e.uid_;
      }
      return graph_ < e.graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for
    // Edge that will not be visible to users, but may be useful within
    // Graph. i.e. Graph needs a way to construct valid Edge objects
    Graph *graph_;         // Pointer back to the Graph
    size_type node1_uid_;  // Node1's uid
    size_type node2_uid_;  // Node2's uid
    size_type uid_;        // This Edge's uid

    /** Private Constructor for a valid edge. */
    Edge(const Graph *graph, const size_type &node1_uid,
         const size_type &node2_uid, const size_type &uid)
        : graph_(const_cast<Graph *>(graph)),
          node1_uid_(node1_uid),
          node2_uid_(node2_uid),
          uid_(uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return i2u_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type uid = i2u_edges_[i];
    return Edge(this, edges_[uid].node1_uid_, edges_[uid].node2_uid_, uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    // HW0: YOUR CODE HERE
    if (!has_node(a) || !has_node(b)) {
      return false;
    }
    // Check if the mapping_ btw a and b exists
    if (mapping_.find(a.uid_) != mapping_.end() &&
        mapping_.at(a.uid_).find(b.uid_) != mapping_.at(a.uid_).end()) {
      // Check that edge is valid
      size_type edge_uid = mapping_.at(a.uid_).at(b.uid_);
      return is_valid_edge(edge_uid);
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already
   * exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges()
   * + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node &a, const Node &b,
                const edge_value_type &value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Add edge if it does not already exist
    if (!has_edge(a, b)) {
      size_type edge_uid = edges_.size();
      size_type edge_idx = i2u_edges_.size();
      i2u_edges_.push_back(edge_uid);
      edges_.push_back(internal_edge(a.uid_, b.uid_, value, edge_idx));
      mapping_[a.uid_][b.uid_] = edge_uid;
      mapping_[b.uid_][a.uid_] = edge_uid;

      return Edge(this, a.uid_, b.uid_, edge_uid);
    }

    // Find edge's index if it already exists
    size_type edge_uid = mapping_[a.uid_][b.uid_];
    return Edge(this, a.uid_, b.uid_, edge_uid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    mapping_.clear();
    edges_.clear();
    nodes_.clear();
    i2u_nodes_.clear();
    i2u_edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    // iterator_category changes to allow NodeIterator support std::min
    // element
    using value_type = Node;                 // Element type
    using pointer = Node *;                  // Pointers to elements
    using reference = Node &;                // Reference to elements
    using difference_type = std::ptrdiff_t;  // Signed difference
    using iterator_category =
        std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference operator for NodeIterator */
    Node operator*() const {
      return Node(graph_, graph_->i2u_nodes_[idx_]);
    }

    /** Increment operator for NodeIterator
     *  Skips a node if it is invalid
     */
    NodeIterator &operator++() {
      ++idx_;
      return *this;
    }

    /** Is equal operator for NodeIterator */
    bool operator==(const NodeIterator &node_iter) const {
      return graph_ == node_iter.graph_ && idx_ == node_iter.idx_;
    }

    /** Not Equal operator for NodeIterator */
    bool operator!=(const NodeIterator &node_iter) const {
      return !((*this) == node_iter);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph *graph_;   // Pointer to the Graph
    size_type idx_;  // Current index

    /** Construct a valid NodeIterator. */
    NodeIterator(const Graph *graph, const size_type &idx)
        : graph_(const_cast<Graph *>(graph)), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns NodeIterator pointing to the first element in the nodes'
   * container. */
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  /** Returns NodeIterator pointing to the past-the-end element in the
   * nodes' container. */
  node_iterator node_end() const { return NodeIterator(this, num_nodes()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge *;                             // Pointers to elements
    using reference = Edge &;                           // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference operator for IncidentIterator */
    Edge operator*() const {
      return Edge(graph_, center_node_, iter_->first, iter_->second);
    }

    /** Increment operator for NodeIterator */
    IncidentIterator &operator++() {
      ++iter_;
      // Skip the edges that were removed
      while (iter_ != end_ && !graph_->is_valid_edge(iter_->second)) {
        ++iter_;
      }
      return *this;
    }

    /** Is equal operator for NodeIterator */
    bool operator==(const IncidentIterator &inc_iter) const {
      return graph_ == inc_iter.graph_ &&
             center_node_ == inc_iter.center_node_ && iter_ == inc_iter.iter_ &&
             end_ == inc_iter.end_;
    }

    /** Is Not equal operator for NodeIterator */
    bool operator!=(const IncidentIterator &inc_iter) const {
      return !((*this) == inc_iter);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    typedef typename std::unordered_map<size_type, size_type>::const_iterator
        map_iterator_;
    Graph *graph_;           // Pointer to the Graph
    size_type center_node_;  // Node that spawns the IncidentIterator (node1)
    // const_iterator for unordered map (node2 -> edge_uid)
    map_iterator_ iter_;
    map_iterator_ end_;

    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Graph *graph, const size_type &center_idx,
                     const map_iterator_ &iter, const map_iterator_ &end)
        : graph_(const_cast<Graph *>(graph)),
          center_node_(center_idx),
          iter_(iter),
          end_(end) {
      // Skip the edges that were removed
      while (iter_ != end_ && !graph_->is_valid_edge(iter_->second)) {
        ++iter_;
      }
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge *;                             // Pointers to elements
    using reference = Edge &;                           // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference operator for EdgeIterator */
    Edge operator*() const {
      size_type edge_uid = graph_->i2u_edges_[idx_];
      return Edge(graph_, graph_->edges_[edge_uid].node1_uid_,
                  graph_->edges_[edge_uid].node2_uid_, edge_uid);
    }

    /** Increment operator for EdgeIterator */
    EdgeIterator &operator++() {
      ++idx_;
      return *this;
    }

    /** Is equal operator for EdgeIterator */
    bool operator==(const EdgeIterator &edge_iter) const {
      return graph_ == edge_iter.graph_ && idx_ == edge_iter.idx_;
    }

    /** Is Not equal operator for EdgeIterator */
    bool operator!=(const EdgeIterator &edge_iter) const {
      return !((*this) == edge_iter);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph *graph_;   // Pointer to the Graph
    size_type idx_;  // Current index

    /** Construct an valid EdgeIterator. */
    EdgeIterator(const Graph *graph, const size_type &idx)
        : graph_(const_cast<Graph *>(graph)), idx_(idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns EdgeIterator at pointing to the first element in the edges'
   * container. */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  /** Returns EdgeIterator referring to the past-the-end element in the
   * edges' container. */
  edge_iterator edge_end() const { return EdgeIterator(this, num_edges()); }

  /** Remove an edge to the graph if it exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if an edge is removed; 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b),
   *             new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   *
   * The internal information won’t be deleted, instead the edge's uid wil be
   * removed from the i2u_eges_ vector, ensuring that the user cannot access
   * that information.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * EdgeIterator and IncidentIterator will skip edges that have been removed.
   *
   * Complexity: O(1).
   */
  size_type remove_edge(const Node &a, const Node &b) {
    // Edge does not exist
    if (!has_edge(a, b)) {
      return 0;
    }
    // Swap the uid of the edge to be removed with the uid of the last element
    // of i2u_edges_ and remove the new last element
    size_type edge_uid = mapping_[a.uid_][b.uid_];
    size_type edge_idx = edges_[edge_uid].idx_;
    i2u_edges_[edge_idx] = i2u_edges_.back();
    i2u_edges_.back() = edge_uid;
    i2u_edges_.pop_back();
    // Update the indicies of swapped edges
    edges_[edge_uid].idx_ = i2u_edges_.size();
    edges_[i2u_edges_[edge_idx]].idx_ = edge_idx;

    return 1;
  }

  /** Remove an edge to the graph if it exists.
   * @pre @a e is an edge object
   * @return 1 if an edge is removed; 0 otherwise
   * @post has_edge(e.node1(), e.node2()) == false
   * @post If old has_edge(e.node1(), e.node2()),
   *             new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   *
   * The internal information won’t be deleted, instead the edge's uid will be
   * removed from the i2u_eges_ vector, ensuring that the user cannot access
   * that information.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * EdgeIterator and IncidentIterator will skip edges that have been removed.
   *
   * Complexity: O(1).
   */
  size_type remove_edge(const Edge &e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge to the graph if it exists.
   * @pre @a e_it is an EdgeIterator object
   * @return EdgeIterator object unchanged
   * @post has_edge((*e_it).node1(), (*e_it).node2()) == false
   * @post If old has_edge((*e_it).node1(), (*e_it).node2()),
   *             new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   *
   * The internal information won’t be deleted, instead the edge's uid will be
   * removed from the i2u_eges_ vector, ensuring that the user cannot access
   * that information.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * EdgeIterator and IncidentIterator will skip edges that have been removed.
   *
   * Complexity: O(1).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove a node to the graph if it exists.
   * @pre @a n_it is a NodeIterator object
   * @return 1 if an node is removed; 0 otherwise
   * @post has_node(*n_it) == false
   * @post If old has_node(*n_it),
   *             new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   *
   * The internal information won’t be deleted, instead the node's uid and uid
   * of its incident edges will be removed from the i2u_nodes_ and i2u_edges_
   * vectors respectively, ensuring that the user cannot access that
   * information.
   *
   * Invalidates node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Node objects.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * NodeIterator will skip nodes that have been removed.
   *
   * EdgeIterator and IncidentIterator will skip edges that have been removed.
   *
   * Complexity: O(1). This assumes that the Graph is sparse, so that the
   * maximum degree of a node is much less than the number of nodes. True
   * complexity is O((*n_it).degree()).
   */
  size_type remove_node(const Node &n) {
    // Node does not exist
    if (!has_node(n)) {
      return 0;
    }
    // Remove incident edges
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it) {
      remove_edge(*it);
    }
    // Swap the node to be removed with the last element of nodes_
    // and remove the new last element - O(1)
    size_type n_idx = nodes_[n.uid_].idx_;
    size_type n_uid = n.uid_;
    i2u_nodes_[n_idx] = i2u_nodes_.back();
    i2u_nodes_.back() = n_uid;
    i2u_nodes_.pop_back();
    // Update indecies of swapped nodes
    nodes_[n_uid].idx_ = i2u_nodes_.size();
    nodes_[i2u_nodes_[n_idx]].idx_ = n_idx;

    return 1;
  }

  /** Remove a node to the graph if it exists.
   * @pre @a n_it is a NodeIterator object
   * @return NodeIterator object unchanged
   * @post has_node(*n_it) == false
   * @post If old has_node(*n_it),
   *             new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   *
   * The internal information won’t be deleted, instead the node's uid and uid
   * of its incident edges will be removed from the i2u_nodes_ and i2u_edges_
   * vectors respectively, ensuring that the user cannot access that
   * information.
   *
   * Invalidates node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Node objects.
   *
   * Invalidates edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * NodeIterator will skip nodes that have been removed.
   *
   * EdgeIterator and IncidentIterator will skip edges that have been removed.
   *
   * Complexity: O(1). This assumes that the Graph is sparse, so that the
   * maximum degree of a node is much less than the number of nodes. True
   * complexity is O((*n_it).degree()).
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // Internal type for Node
  struct internal_node {
    Point position_;         //< Node's position
    node_value_type value_;  //< Node's value
    size_type idx_;          //< Node's index in i2u_nodes_

    internal_node(const Point &pos, const node_value_type &val,
                  const size_type &idx)
        : position_(pos), value_(val), idx_(idx) {}
  };

  // Internal type for Edge
  struct internal_edge {
    size_type node1_uid_;    //< Node1's uid
    size_type node2_uid_;    //< Node2's uid
    edge_value_type value_;  //< Edge's value
    size_type idx_;          //< Edge's index in i2u_edges_

    internal_edge(const size_type &uid1, const size_type &uid2,
                  const edge_value_type &val, const size_type &idx)
        : node1_uid_(uid1), node2_uid_(uid2), value_(val), idx_(idx) {}
  };

  /* Stores internal_node for any node added, even if later removed.
   * Indexed by node uid.  */
  std::vector<internal_node> nodes_;
  /* Store the currently "active" set of nodes.
   * Indexed by node idx.  */
  std::vector<size_type> i2u_nodes_;

  /* Stores internal_edge for any edge edges, even if later removed.
   * Indexed by edge uid.  */
  std::vector<internal_edge> edges_;
  /* Store the currently "active" set of nodes.
   * Indexed by edge idx.  */
  std::vector<size_type> i2u_edges_;

  /* Maps node1_uid -> node2_uid -> edge_uid */
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>>
      mapping_;

  /** Test whether Node is valid based on its uid. */
  bool is_valid_node(const size_type &uid) const {
    return uid >= 0 && uid < nodes_.size() &&
           nodes_[uid].idx_ < i2u_nodes_.size();
  }

  /** Test whether Edge is valid based on its uid. */
  bool is_valid_edge(const size_type &uid) const {
    return uid >= 0 && uid < edges_.size() &&
           edges_[uid].idx_ < i2u_edges_.size();
  }
};

#endif  // CME212_GRAPH_HPP
