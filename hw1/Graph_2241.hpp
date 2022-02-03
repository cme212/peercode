#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
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
template <typename V>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Synonym for V. */
  using node_value_type = V;

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
      return graph_->nodes_[idx_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value, of type V. */
    node_value_type &value() { return graph_->nodes_[idx_].value_; }

    /** Return this node's value, of type V; const version */
    const node_value_type &value() const { return graph_->nodes_[idx_].value_; }

    /** Return number of edges connected to this node */
    size_type degree() const {
      // Check if a Node has edges incident to it
      bool connected = graph_->mapping_.find(idx_) != graph_->mapping_.end();
      return connected ? graph_->mapping_.at(idx_).size() : 0;
    }

    /** Returns IncidentIterator at pointing to the first element of
     * incident edges to this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, idx_, graph_->mapping_[idx_].begin());
    }
    /** Returns IncidentIterator referring to the past-the-end element
     * of incident edges to this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, idx_, graph_->mapping_[idx_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (idx_ == n.idx_);
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
        return idx_ < n.idx_;
      }
      return graph_ < n.graph_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for
    // Node that will not be visible to users, but may be useful within
    // Graph. i.e. Graph needs a way to construct valid Node objects
    Graph *graph_;   // Pointer back to the Graph
    size_type idx_;  // This Node's index

    /** Private Constructor for a valid node. */
    Node(const Graph *graph, const size_type &idx)
        : graph_(const_cast<Graph *>(graph)), idx_(idx) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return nodes_.size(); }

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
    nodes_.push_back(internal_node(position, value, true));
    return Node(this, nodes_.size() - 1);  // Valid Node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_ && n.idx_ < num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);  // Valid Node
  }

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
    Edge() : graph_(nullptr), valid_(0) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      // HW0: YOUR CODE HERE
      bool exact_order =
          node1_idx_ == e.node1_idx_ && node2_idx_ == e.node2_idx_;
      bool flipped_order =
          node1_idx_ == e.node2_idx_ && node2_idx_ == e.node1_idx_;
      return graph_ == e.graph_ && (exact_order || flipped_order);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // HW0: YOUR CODE HERE
      if (graph_ == e.graph_) {
        return node1_idx_ < e.node1_idx_ && node2_idx_ < e.node2_idx_;
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
    size_type node1_idx_;  // Index of node1
    size_type node2_idx_;  // Index of node2
    bool valid_;           // Edge's valid status

    /** Private Constructor for a valid edge. */
    Edge(const Graph *graph, const size_type &node1_idx,
         const size_type &node2_idx, const bool &valid)
        : graph_(const_cast<Graph *>(graph)),
          node1_idx_(node1_idx),
          node2_idx_(node2_idx),
          valid_(valid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    // HW0: YOUR CODE HERE
    return mapping_.find(a.idx_) != mapping_.end() &&
           mapping_.at(a.idx_).find(b.idx_) != mapping_.at(a.idx_).end();
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
  Edge add_edge(const Node &a, const Node &b) {
    // HW0: YOUR CODE HERE
    Edge new_edge = Edge(this, a.idx_, b.idx_, 1);
    // add edge if it does not already exist
    if (!has_edge(a, b)) {
      edges_.push_back(new_edge);
      size_type edge_idx = edges_.size() - 1;
      mapping_[a.idx_][b.idx_] = edge_idx;
      mapping_[b.idx_][a.idx_] = edge_idx;
    }
    return new_edge;
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
    Node operator*() const { return Node(graph_, node_idx_); }

    /** Increment operator for NodeIterator */
    NodeIterator &operator++() {
      ++node_idx_;
      return *this;
    }

    /** Is equal operator for NodeIterator */
    bool operator==(const NodeIterator &node_iter) const {
      return graph_ == node_iter.graph_ && node_idx_ == node_iter.node_idx_;
    }

    /** Not Equal operator for NodeIterator */
    bool operator!=(const NodeIterator &node_iter) const {
      return !((*this) == node_iter);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph *graph_;        // Pointer to the Graph
    size_type node_idx_;  // Position in the node's vector

    /** Construct a valid NodeIterator. */
    NodeIterator(const Graph *graph, const size_type &idx)
        : graph_(const_cast<Graph *>(graph)), node_idx_(idx) {}
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
      return Edge(graph_, center_node_, iter_->first, 1);
    }

    /** Increment operator for NodeIterator */
    IncidentIterator &operator++() {
      ++iter_;
      return *this;
    }

    /** Is equal operator for NodeIterator */
    bool operator==(const IncidentIterator &inc_iter) const {
      return graph_ == inc_iter.graph_ &&
             center_node_ == inc_iter.center_node_ && iter_ == inc_iter.iter_;
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
    // const_iterator for unordered map (node2 -> edge_idx)
    map_iterator_ iter_;

    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Graph *graph, const size_type &center_idx,
                     const map_iterator_ &iter)
        : graph_(const_cast<Graph *>(graph)),
          center_node_(center_idx),
          iter_(iter) {}
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
    Edge operator*() const { return graph_->edges_[edge_idx_]; }

    /** Increment operator for EdgeIterator */
    EdgeIterator &operator++() {
      ++edge_idx_;
      return *this;
    }

    /** Is equal operator for EdgeIterator */
    bool operator==(const EdgeIterator &edge_iter) const {
      return graph_ == edge_iter.graph_ && edge_idx_ == edge_iter.edge_idx_;
    }

    /** Is Not equal operator for EdgeIterator */
    bool operator!=(const EdgeIterator &edge_iter) const {
      return !((*this) == edge_iter);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph *graph_;        // Pointer to the Graph
    size_type edge_idx_;  // Position in the edges's vector

    /** Construct an valid EdgeIterator. */
    EdgeIterator(const Graph *graph, const size_type &idx)
        : graph_(const_cast<Graph *>(graph)), edge_idx_(idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns EdgeIterator at pointing to the first element in the edges'
   * container. */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  /** Returns EdgeIterator referring to the past-the-end element in the
   * edges' container. */
  edge_iterator edge_end() const { return EdgeIterator(this, num_edges()); }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // Internal type for Node
  struct internal_node {
    Point position_;
    node_value_type value_;
    bool valid_;

    internal_node(const Point &pos, const node_value_type &val,
                  const bool &valid)
        : position_(pos), value_(val), valid_(valid) {}

    internal_node()
        : position_(Point()), value_(node_value_type()), valid_(false) {}
  };

  // This Graph's Nodes
  std::vector<internal_node> nodes_;
  // This Graph's Edges
  std::vector<Edge> edges_;
  // Mapping of nodes' connectedness to edges' unique id
  // node1_idx -> node2_idx -> edge_idx
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>>
      mapping_;
};

#endif  // CME212_GRAPH_HPP
