#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V> class Graph {
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

  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph *graph_;
  Graph() {
    points_ = {};
    nodes_ = {};
    edges_ = {};
    values_ = {};
    edgeSet_ = {};
    incidents_ = {};
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
    // Graph *graph_;
    Node() {
      idx_ = 0;
      graph_ = NULL;
      graph_->values_[idx_] = node_value_type();
    }

    /** Return this node's position. */
    const Point &position() const { return graph_->points_[idx_]; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return idx_; }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type &value() { return graph_->values_[idx_]; }

    const node_value_type &value() const { return graph_->values_[idx_]; }

    size_type degree() const { return (graph_->incidents_[idx_]).size(); }

    /** begin() function for incident iterator. */
    incident_iterator edge_begin() const {
      incident_iterator new_it;
      typename std::vector<Edge>::const_iterator it =
          (graph_->incidents_[idx_]).begin();
      new_it.curIncident = it;
      return new_it;
    }

    /** end() function for incident iterator. */
    incident_iterator edge_end() const {
      incident_iterator new_it;
      typename std::vector<Edge>::const_iterator it =
          (graph_->incidents_[idx_]).end();
      new_it.curIncident = it;
      return new_it;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      if (graph_ == n.graph_ && idx_ == n.idx_) {
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
    bool operator<(const Node &n) const { return idx_ < n.idx_; }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type idx_;
    Graph *graph_;

    Node(size_type n, Graph *graph) {
      idx_ = n;
      graph_ = graph;
      graph_->values_[idx_] = node_value_type();
    }
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
   * @param[in] v The new node's value (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type& v = node_value_type()) {
    Node new_node = Node(num_nodes(), this);
    assert(sizeof(new_node) <= 16);
    points_.push_back(position);
    nodes_.push_back(new_node);
    incidents_[nodes_.size()] = {};
    values_[nodes_.size()] = v;
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    if (this == n.graph_ && n.idx_ < size()) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    return nodes_[i];
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
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() {
      n1_ = NULL;
      n2_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const { return *n1_; }

    /** Return the other node of this Edge */
    Node node2() const { return *n2_; }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      return (((n1_->idx_ == e.n1_->idx_) && (n2_->idx_ == e.n2_->idx_)) ||
              ((n1_->idx_ == e.n2_->idx_) && (n2_->idx_ == e.n1_->idx_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      size_type min1 = std::min(n1_->idx_, n2_->idx_);
      size_type max1 = std::max(n1_->idx_, n2_->idx_);
      size_type min2 = std::min(e.n1_->idx_, e.n2_->idx_);
      size_type max2 = std::max(e.n1_->idx_, e.n2_->idx_);
      if (min1 < min2) {
        return true;
      } else if (min1 > min2) {
        return false;
      } else {
        return (max1 < max2);
      }
    }

    /** Return the other node of the edge */
    Node otherNode(Node given_node) {
      Node x = *n1_, y = *n2_;

      return x == given_node ? y : x;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const Node *n1_, *n2_;
    Edge(const Node &a, const Node &b) {
      n1_ = &a;
      n2_ = &b;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return edges_.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    size_type x = std::min(a.idx_, b.idx_), y = std::max(a.idx_, b.idx_);
    std::string edge_hash = std::to_string(x) + '-' + std::to_string(y);
    return (edgeSet_.find(edge_hash) != edgeSet_.end());
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
    Edge new_edge(a, b);
    assert(sizeof(new_edge) <= 32);
    if (!has_edge(a, b)) {
      size_type x = std::min(a.idx_, b.idx_), y = std::max(a.idx_, b.idx_);
      edges_.push_back(new_edge);
      // create a string hash for edgeSet_
      edgeSet_.insert(std::to_string(x) + '-' + std::to_string(y));
      incidents_[a.idx_].push_back(new_edge);
      if (a.idx_ != b.idx_) { // In case the nodes are the same
        Edge flipped_edge(b, a);
        incidents_[b.idx_].push_back(flipped_edge);
      }
    }
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points_ = {};
    nodes_ = {};
    edges_ = {};
    edgeSet_ = {};
    incidents_ = {};
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    Node operator*() const { return *curNode; };
    NodeIterator &operator++() {
      NodeIterator next_node_it;
      next_node_it.curNode = curNode++;
      return next_node_it;
    };
    bool operator==(const NodeIterator &nodeIter) const {
      return curNode == nodeIter.curNode;
    };

  private:
    friend class Graph;
    NodeIterator(typename std::vector<Node>::const_iterator it) {
      curNode = it;
    }
    typename std::vector<Node>::const_iterator curNode;
  };

  node_iterator node_begin() const {
    node_iterator new_it;
    typename std::vector<Node>::const_iterator it = nodes_.begin();
    new_it.curNode = it;
    return new_it;
  };
  node_iterator node_end() const {
    node_iterator new_it;
    typename std::vector<Node>::const_iterator it = nodes_.end();
    new_it.curNode = it;
    return new_it;
  };

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    Edge operator*() const { return *curIncident; };

    IncidentIterator &operator++() {
      IncidentIterator it;
      it.curIncident = curIncident++;
      return it;
    };

    bool operator==(const IncidentIterator &IncidentIter) const {
      return curIncident == IncidentIter.curIncident;
    };

  private:
    friend class Graph;
    typename std::vector<Edge>::const_iterator curIncident;
    IncidentIterator(typename std::vector<Edge>::const_iterator it) {
      curIncident = it;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    Edge operator*() const { return *curEdge; };
    EdgeIterator &operator++() {
      EdgeIterator it;
      it.curEdge = curEdge++;
      return it;
    };
    bool operator==(const EdgeIterator &EdgeIter) const {
      return curEdge == EdgeIter.curEdge;
    };

  private:
    friend class Graph;
    EdgeIterator(typename std::vector<Edge>::const_iterator it) {
      curEdge = it;
    }
    typename std::vector<Edge>::const_iterator curEdge;
  };

  /** begin() function for edge iterator. */
  edge_iterator edge_begin() const {
    edge_iterator new_it;
    typename std::vector<Edge>::const_iterator it = edges_.begin();
    new_it.curEdge = it;
    return new_it;
  };
  
  /** end() function for edge iterator. */
  edge_iterator edge_end() const {
    edge_iterator new_it;
    typename std::vector<Edge>::const_iterator it = edges_.end();
    new_it.curEdge = it;
    return new_it;
  };

private:
  std::vector<Point> points_;
  std::vector<Node> nodes_;
  std::vector<Edge> edges_;
  std::unordered_map<int, node_value_type> values_;
  std::unordered_set<std::string> edgeSet_;
  std::unordered_map<size_type, std::vector<Edge>> incidents_;
};

#endif // CME212_GRAPH_HPP
