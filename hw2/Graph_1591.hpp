#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int>
class Graph {
private:
  // Forward declaration of internal type for nodes
  struct internal_node;

public:
  // ==========================================================================
  // PUBLIC TYPE DEFINITIONS
  // ==========================================================================

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Node value type */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Edge value type */
  using edge_value_type = E;

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

  // ==========================================================================
  // CONSTRUCTORS AND DESTRUCTOR
  // ==========================================================================

  /** Construct an empty graph. */
  Graph() : nodes_(), valid_uids_(), num_edges_(0), adj_uids_max_(), adj_uids_min_(), edge_values_() {};

  /** Default destructor */
  ~Graph() = default;

  // ==========================================================================
  // NODES
  // ==========================================================================

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
  public:
    /** Construct an invalid Node.
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
    Node() : graph_(nullptr), uid_(0) {}

    /** Retrieve the position of this node.
     *
     * @return  A reference to a Point holding this node's position.
     *
     * @post    The Point cannot be modified via this reference.
     *
     * The complexity of retrieval is O(1).
     */
    const Point& position() const {
      assert(this->graph_ != nullptr);  // this node must belong to a graph
      return graph_->nodes_[uid_].point;
    }

    /** Retrieve the position of this node.
     *
     * @return  A reference which can be used to modify this node's
     *          position, i.e., Point.
     *
     * The complexity of retrieval is O(1).
     */
    Point& position() {
      return const_cast<Point&>(const_cast<const Node*>(this)->position());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(this->graph_ != nullptr);  // this node must belong to a graph
      return graph_->nodes_[uid_].index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:


    /** Retrieve the value of this node.
     *
     * @return  A reference to this node's internal value.
     *
     * @post    The node's internal value cannot be modified via this reference.
     *
     * The complexity of retrieval is O(1).
     */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value;
    }

    /** Retrieve the value of this node.
     *
     * @return  A reference which can be used to modify this node's internal value.
     *
     * The complexity of retrieval is O(1).
     */
    node_value_type& value() {
      return const_cast<node_value_type&>(const_cast<const Node*>(this)->value());
    }

    /** Compute the number of edges incident to this node.
     *
     * @return  The number of edges incident to this node.
     *
     * @post    0 <= degree() < graph_->num_nodes().
     *
     * The complexity of this operation is O(1).
     */
    size_type degree() const {
      // Add up the number of "max edges" and "min edges" incident to this node
      return graph_->adj_uids_max_[uid_].size() + graph_->adj_uids_min_[uid_].size();
    }

    /** Return an iterator pointed at the first edge incident to this node.
     *
     * @return  An IncidentIterator pointed to the first edge incident to this node.
     *
     * @pre     The node spawning this iterator is valid
     *          (i.e., graph_ != nullptr and there exists 0 <= i < graph_->num_nodes()
     *           such that uid_ == graph_->valid_uids_[i]).
     *
     * The complexity of this operation is O(1). If there are no edges incident
     * to this node, the returned iterator equals edge_end().
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /** Return an iterator pointed at one past the last edge incident to this node.
     *
     * @return  An IncidentIterator pointed at one past the last edge incident to this node.
     *          Attempting to dereference this iterator results in undefined behaviour.
     *
     * @pre     The node spawning this iterator is valid
     *          (i.e., graph_ != nullptr and there exists 0 <= i < graph_->num_nodes()
     *           such that uid_ == graph_->valid_uids_[i]).
     *
     * The complexity of this operation is O(1). If there are no edges incident
     * to this node, the returned iterator equals edge_begin().
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Due to how uids are created (see add_node), we can compare them in
      // place of indices
      return graph_ == n.graph_ && uid_ == n.uid_;
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
      return std::tie(graph_, uid_) < std::tie(n.graph_, n.uid_);
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer to this node's graph (can be nullptr if this node is invalid)
    graph_type* graph_;

    // Unique identifier used to access this node's data in the graph
    size_type uid_;

    /** Construct a valid Node.
     *
     *  @param[in] graph  The graph containing this node.
     *  @param[in] uid    The unique identifier for this node in the graph.
     *
     *  The complexity of this constructor is O(1).
     */
    Node(const Graph* graph, size_type uid) :
      graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  // END NODE =================================================================

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return valid_uids_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @pre  node_value_type has a valid default constructor.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    nodes_.emplace_back(position, value, valid_uids_.size());
    valid_uids_.push_back(nodes_.size() - 1);
    adj_uids_min_.resize(nodes_.size());
    adj_uids_max_.resize(nodes_.size());
    edge_values_.resize(nodes_.size());
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this                                // correct graph pointer
           && n.index() < n.graph_->size()                 // node index is in-bounds
           && n.uid_ == n.graph_->valid_uids_[n.index()];  // node uid matches the internal uid at the same index
  }

  /**
   * @brief Return the node with index @a i.
   *
   * @param[in] i   Index of the node to return.
   *
   * @pre   0 <= @a i < num_nodes()
   * @post  result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    return Node(this, valid_uids_[i]);
  }

  // ==========================================================================
  // EDGES
  // ==========================================================================

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr), uid1_(0), uid2_(0) {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_ == e.graph_ && ((uid1_ == e.uid1_ && uid2_ == e.uid2_)
                                    || (uid1_ == e.uid2_ && uid2_ == e.uid1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // We need to be careful of the case when comparing two edges with the
      // same nodes in opposite order; to preserve trichotomy, we first order
      // the nodes for an edge according to their uids before comparing
      size_type uid_min = std::min(uid1_, uid2_);
      size_type uid_max = std::max(uid1_, uid2_);
      size_type e_uid_min = std::min(e.uid1_, e.uid2_);
      size_type e_uid_max = std::max(e.uid1_, e.uid2_);
      return std::tie(graph_, uid_min, uid_max) < std::tie(e.graph_, e_uid_min, e_uid_max);
    }

    /** Compute the length of this Edge.
     *
     * @return  The length of this edge as the Euclidean distance between its
     *          nodes' positions.
     *
     * @pre     The nodes of this edge have valid positions described as Points.
     *
     * The complexity of this operation is O(1).
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Retrieve the value of this edge.
     *
     * @return  A reference to this edge's value.
     *
     * @post    The edge's value cannot be modified via this reference.
     *
     * The worst-case complexity of retrieval is O(max(node1().degree(), node2().degree())).
     */
    const edge_value_type& value() const {
      // Get the unique IDs of the "min" and "max" nodes
      size_type uid_min = std::min(uid1_, uid2_);
      size_type uid_max = std::max(uid1_, uid2_);

      // Find the "max edge" index into edge_values_
      size_type idx = 0;
      for (std::vector<size_type>::const_iterator it = graph_->adj_uids_max_[uid_min].begin();
           it != graph_->adj_uids_max_[uid_min].end();
           ++it) {
        if (*it == uid_max) break;
        ++idx;
      }
      return graph_->edge_values_[uid_min][idx];
    }

    /** Retrieve the value of this edge.
     *
     * @return  A reference which can be used to modify this edge's value.
     *
     * The worst-case complexity of retrieval is O(max(node1().degree(), node2().degree())).
     */
    edge_value_type& value() {
      return const_cast<edge_value_type&>(const_cast<const Edge*>(this)->value());
    }


  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer to this edge's graph (can be nullptr if this edge is invalid)
    graph_type* graph_;

    // Unique identifiers used to access the data for each node of this edge
    // in the graph
    size_type uid1_;
    size_type uid2_;

    /** Construct a valid Edge.
     *  @param[in] graph  The graph containing this edge.
     *  @param[in] uid1   The unique identifier for the first node in the graph.
     *  @param[in] uid2   The unique identifier for the second node in the graph.
     *  Complexity: O(1).
     */
    Edge(const Graph* graph, size_type uid1, size_type uid2) :
      graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {}
  };

  // END EDGE =================================================================

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());  // avoid out-of-bounds issues
    size_type uid1;
    size_type adj_index = i;
    for (uid1 = 0; uid1 < adj_uids_max_.size(); ++uid1) {
      // If adj_index is too large to index into the vector of nodes adjacent
      // to the current node, decrement adj_index by the number of adjacencies
      // and continue to the next node
      if (adj_index >= adj_uids_max_[uid1].size()) {
        adj_index -= adj_uids_max_[uid1].size();
      } else {
      // Otherwise, we have found the node and the index of its adjacent node
      // that together form the i-th edge
        break;
      }
    }
    return Edge(this, uid1, adj_uids_max_[uid1][adj_index]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(max(_a_.degree(), _b_.degree())).
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a == b) {
      return false;  // no edges from a node to itself
    }
    size_type uid_min = std::min(a.uid_, b.uid_);
    size_type uid_max = std::max(a.uid_, b.uid_);
    for (size_type i = 0; i < adj_uids_max_[uid_min].size(); ++i) {
      if (adj_uids_max_[uid_min][i] == uid_max) return true;
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
    assert(!(a == b));  // no edges can be added from a node to itself
    if (!has_edge(a, b)) {
      size_type uid_min = std::min(a.uid_, b.uid_);
      size_type uid_max = std::max(a.uid_, b.uid_);
      adj_uids_max_[uid_min].push_back(uid_max);
      adj_uids_min_[uid_max].push_back(uid_min);
      edge_values_[uid_min].push_back(edge_value_type());
      ++num_edges_;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    valid_uids_.clear();
    num_edges_ = 0;
    adj_uids_max_.clear();
    adj_uids_min_.clear();
    edge_values_.clear();
  }

  // ==========================================================================
  // NODE ITERATOR
  // ==========================================================================

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), idx_(0) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Access the node at this iterator's current position via dereferencing.
     *
     * @return  A Node object corresponding to this iterator's current position.
     *
     * @pre     This iterator != node_end() of the underlying graph.
     * @post    The returned Node cannot be modified.
     *
     * The complexity of this operation is O(1).
     */
    Node operator*() const {
      return Node(graph_, graph_->valid_uids_[idx_]);
    }

    /** Pre-increment this iterator to the next node.
     *
     * @return  A reference to this iterator after it has been incremented.
     *
     * @pre     This iterator != node_end() of the underlying graph.
     * @post    The @a idx_ of this iterator is incremented by one.
     *
     * The complexity of this operation is O(1).
     */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Compare this iterator to another NodeIterator.
     *
     * @param[in] it  A reference to the other NodeIterator to be compared.
     * @return        True if _it_ has the same graph pointer and node index
     *                as this iterator, otherwise False.
     *
     * @pre     Each iterator was spawned by a valid graph
     *          (i.e., graph_ and it.graph_ are not null pointers).
     * @pre     0 <= idx_ < graph_->num_nodes() and 0 <= it.idx_ < it.graph_->num_nodes().
     * @post    Neither iterator is modified.
     *
     * The complexity of this operation is O(1).
     */
    bool operator==(const NodeIterator& it) const {
      return graph_ == it.graph_ && idx_ == it.idx_;
    }


  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;   // pointer to the current node's graph
    size_type idx_;  // index of the current node

    NodeIterator(const Graph* graph, size_type idx) :
      graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator pointed at the first node in this graph.
   *
   * @return  A NodeIterator pointed at the first node in this graph.
   *
   *
   * The complexity of this operation is O(1). If there are no nodes in this
   * graph, the returned iterator equals node_end().
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an iterator pointed at one past the last node in this graph.
   *
   * @return  A NodeIterator pointed at one past the last node in this graph.
   *          Attempting to dereference this iterator results in undefined behaviour.
   *
   * The complexity of this operation is O(1). If there are no edges incident
   * to this node, the returned iterator equals node_begin().
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  // ==========================================================================
  // INCIDENT ITERATOR
  // ==========================================================================

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr), uid1_(0), idx_incident_(0) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Access the incident edge at this iterator's current position via dereferencing.
     *
     * @return  An Edge object corresponding to this iterator's current position.
     *
     * @pre     This iterator != edge_end() of the underlying node.
     * @post    The returned Edge cannot be modified.
     * @post    The node that spawned this iterator == edge.node1().
     *
     * The complexity of this operation is O(1).
     */
    Edge operator*() const {
      size_type uid2;
      size_type num_min_edges = graph_->adj_uids_min_[uid1_].size();
      if (idx_incident_ < num_min_edges) {
        uid2 = graph_->adj_uids_min_[uid1_][idx_incident_];
      } else {
        uid2 = graph_->adj_uids_max_[uid1_][idx_incident_ - num_min_edges];
      }
      return Edge(graph_, uid1_, uid2);
    }

    /** Pre-increment this iterator to the next incident edge.
     *
     * @return  A reference to this iterator after it has been incremented.
     *
     * @pre     This iterator != edge_end() of the underlying node.
     * @post    The @a idx_incident_ of this iterator is incremented by one.
     *
     * The complexity of this operation is O(1).
     */
    IncidentIterator& operator++() {
      ++idx_incident_;
      return *this;
    }

    /** Compare this iterator to another IncidentIterator.
     *
     * @param[in] it  A reference to the other IncidentIterator to be compared.
     * @return        True if _it_ has the same graph pointer, node unique id,
     *                and incident index as this iterator, otherwise False.
     *
     * @pre     Each iterator was spawned by a valid node
     *          (i.e., for each iterator, graph_ != nullptr and there exists
     *           0 <= i < graph_->num_nodes() such that uid1_ == graph_->valid_uids_[i]).
     * @post    Neither iterator is modified.
     *
     * The complexity of this operation is O(1).
     */
    bool operator==(const IncidentIterator& it) const {
      return graph_ == it.graph_ && uid1_ == it.uid1_ && idx_incident_ == it.idx_incident_;
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;            // pointer to this iterator's graph
    size_type uid1_;          // the unique id of the node ("node1" of any incident edge)
    size_type idx_incident_;  // the index for incident edges

    IncidentIterator(const Graph* graph, size_type uid1, size_type idx_incident) :
        graph_(const_cast<Graph*>(graph)), uid1_(uid1), idx_incident_(idx_incident) {}
  };

  // ==========================================================================
  // EDGE ITERATOR
  // ==========================================================================

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), idx1_(0), idx_adj_(0) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Access the edge at this iterator's current position via dereferencing.
     *
     * @return  An Edge object corresponding to this iterator's current position.
     *
     * @pre     This iterator != edge_end() of the underlying graph.
     * @post    The returned Edge cannot be modified.
     * @post    edge.node1().index() < edge.node2().index().
     *
     * The complexity of this operation is O(1).
     */
    Edge operator*() const {
      size_type uid1 = graph_->valid_uids_[idx1_];
      size_type uid2 = graph_->adj_uids_max_[uid1][idx_adj_];
      return Edge(graph_, uid1, uid2);
    }

    /** Pre-increment this iterator to the next edge in the graph.
     *
     * @return  A reference to this iterator after it has been incremented.
     *
     * @pre     This iterator != edge_end() of the underlying graph.
     * @post    If this iterator != edge_end(), then @a idx1_ of
     *          this iterator == the index of a valid node degree() >= 1.
     */
    EdgeIterator& operator++() {
      assert(idx1_ < graph_->num_nodes());
      size_type uid1 = graph_->valid_uids_[idx1_];
      assert(!graph_->adj_uids_max_[uid1].empty());
      if (idx_adj_ < graph_->adj_uids_max_[uid1].size() - 1) {
        ++idx_adj_;
      } else {
        // Find the next node with incident edges; if there are none, idx1_
        // reaches graph_->num_nodes(), which makes this iterator equivalent to
        // the end iterator
        idx_adj_ = 0;
        ++idx1_;
        for ( ; idx1_ < graph_->num_nodes(); ++idx1_) {
          uid1 = graph_->valid_uids_[idx1_];
          if (!graph_->adj_uids_max_[uid1].empty()) { break; }
        }
      }
      return *this;
    }

    /** Compare this iterator to another EdgeIterator.
     *
     * @param[in] it  A reference to the other EdgeIterator to be compared.
     * @return        True if _it_ has the same graph pointer, node1 index,
     *                and adjacent node index as this iterator, otherwise False.
     *
     * @pre     Each iterator was spawned by a valid graph
     *          (i.e., for each iterator, graph_ != nullptr and 0 <= idx1_ < graph_->num_nodes()).
     * @post    Neither iterator is modified.
     *
     * The complexity of this operation is O(1).
     */
    bool operator==(const EdgeIterator& it) const {
      return graph_ == it.graph_ && idx1_ == it.idx1_ && idx_adj_ == it.idx_adj_;
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;       // pointer to this iterator's graph
    size_type idx1_;     // the index of the edge node with the smaller unique ID
    size_type idx_adj_;  // index into adjacent nodes along "max edges"

    EdgeIterator(const Graph* graph, size_type idx1, size_type idx_adj) :
        graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx_adj_(idx_adj) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return an iterator pointed at the first edge of this graph.
   *
   * @return  An EdgeIterator pointed to the first edge of this graph.
   *
   * If there are no edges in this graph (i.e., num_edges() == 0), the returned
   * iterator equals edge_end().
   */
  edge_iterator edge_begin() const {
    // Find the first node with a connected "max edge"
    size_type idx1;
    size_type uid1;
    for (idx1 = 0; idx1 < num_nodes(); ++idx1) {
      uid1 = valid_uids_[idx1];
      if (!adj_uids_max_[uid1].empty()) { break; }
    }
    return EdgeIterator(this, idx1, 0);
  }

  /** Return an iterator pointed at one past the last edge of this graph.
   *
   * @return  An EdgeIterator pointed to one past the last edge of this graph.
   *          Attempting to dereference this iterator results in undefined behaviour.
   *
   * The complexity of this operation is O(1). If there are no edges in this
   * graph (i.e., num_edges() == 0), the returned iterator equals edge_begin().
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }

  // ==========================================================================
  // REMOVAL METHODS
  // ==========================================================================

  /** Remove a node and all of its incident edges from this graph.
   *
   * @param[in] n   The Node to remove from this graph.
   *
   * @return    0 if _n_ is not a node in this graph, otherwise 1 to indicate
   *            successful node removal.
   *
   * @post  has_node(_n_) == false.
   * @post  If has_node(_n_) beforehand, then num_nodes() is reduced by one.
   * @post  If has_node(_n_) beforehand, then num_edges() is reduced by _n_.degree().
   * @post  If !has_node(_n_) beforehand, then num_nodes() remains the same.
   * @post  If !has_node(_n_) beforehand, then num_edges() remains the same.
   * @post  If has_node(_n_) beforehand, then any Edge e with e.node1() == _n_
   *        or e.node2() == _n_ has been removed from the graph and is now invalid.
   * @post  All pre-existing NodeIterators are invalidated by this operation.
   * @post  All pre-existing EdgeIterators are invalidated by this operation.
   * @post  Any pre-existing IncidentIterators associated with _n_ are invalidated
   *        by this operation.
   *
   * The worst-case complexity of this operation is O(n.degree()^2), which is
   * much smaller than O(num_nodes()) assuming the graph is sparse.
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
        return 0;
    } else {
        size_type uid = n.uid_;
        size_type idx = n.index();

        // Clear all connected edges
        // Complexity: O(n.degree()^2)
        while (!adj_uids_max_[uid].empty()) {
          remove_edge(n, Node(this, adj_uids_max_[uid][0]));
        }
        while (!adj_uids_min_[uid].empty()) {
          remove_edge(n, Node(this, adj_uids_min_[uid][0]));
        }

        // Update valid unique IDs and internal node data by:
        // 1) moving the last unique ID to the current index
        // 2) enforcing the node index / unique invariant
        // Complexity: O(1)
        valid_uids_[idx] = std::move(valid_uids_.back());
        valid_uids_.pop_back();
        nodes_[valid_uids_[idx]].index = idx;

        return 1;
    }
  }


  /** Remove a node and all of its incident edges from this graph.
   *
   * @param[in] n_it   A NodeIterator corresponding to the Node to remove.
   *
   * @return    node_begin(), valid for this graph after node removal.
   *
   * @pre   n_it is a valid NodeIterator for this graph.
   * @pre   n_it != node_end().
   * @post  num_nodes() is reduced by one.
   * @post  num_edges() is reduced by (*n_it).degree() (dereferenced before removal).
   * @post  Any Edge e with e.node1() == *n_it or e.node2() == *n_it has been
   *        removed from the graph and is now invalid.
   * @post  All pre-existing NodeIterators are invalidated by this operation.
   * @post  All pre-existing EdgeIterators are invalidated by this operation.
   * @post  Any pre-existing IncidentIterators associated with *n_it are
   *        invalidated by this operation.
   *
   * The worst-case complexity of this operation is O(num_nodes()).
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
  }


  /** Remove an edge from this graph.
   *
   * @param[in] a   A Node of the edge to remove.
   * @param[in] b   The other Node of the edge to remove.
   *
   * @return    0 if has_edge(a, b), otherwise 1 to indicate successful edge removal.
   *
   * @pre   has_node(_a_) == true and has_node(_b_) == true.
   * @post  If has_edge(_a_, _b_) beforehand, then num_edges() is reduced by one.
   * @post  If !has_edge(_a_, _b_) beforehand, then num_edges() remains the same.
   * @post  All pre-existing EdgeIterators are invalidated by this operation.
   * @post  All pre-existing IncidentIterators are invalidated by this operation.
   *
   * The worst-case complexity of this operation is O(_a_.degree() + _b_.degree()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) return 0;

    size_type uid_min = std::min(a.uid_, b.uid_);
    size_type uid_max = std::max(a.uid_, b.uid_);

    // Find the adjacency index for the "max" edge orientation
    bool has_edge_ab = false;
    size_type idx_max = 0;
    for ( ; idx_max < adj_uids_max_[uid_min].size(); ++idx_max) {
      if (adj_uids_max_[uid_min][idx_max] == uid_max) {
        has_edge_ab = true;
        break;
      }
    }

    // If we found a matching unique ID, this edge exists and must be deleted
    if (has_edge_ab) {
      // Find the adjacency index for the "min" edge orientation, which has to
      // exist if we already found the "max" index
      size_type idx_min = 0;
      for ( ; idx_min < adj_uids_min_[uid_max].size(); ++idx_min) {
        if (adj_uids_min_[uid_max][idx_min] == uid_min) {
          break;
        }
      }
      assert(idx_min < adj_uids_min_[uid_max].size());

      // Erase the entries corresponding to this edge from the adjacency vectors
      // and edge values vector, thereby moving subsequent entries up by one element
      adj_uids_min_[uid_max].erase(adj_uids_min_[uid_max].begin() + idx_min);
      adj_uids_max_[uid_min].erase(adj_uids_max_[uid_min].begin() + idx_max);
      edge_values_[uid_min].erase(edge_values_[uid_min].begin() + idx_max);

      // Finish off by reducing the number of edges in the graph
      --num_edges_;
    }

    return 1;
  }


  /** Remove an edge from this graph.
   *
   * @param[in] e   The Edge to remove.
   *
   * @return    0 if has_edge(a, b), otherwise 1 to indicate successful edge removal.
   *
   * @pre   has_node(_e_.node1()) == true and has_node(_e_.node2()) == true.
   * @post  If has_edge(_e_.node1(), _e_.node2()), then num_edges() is reduced by one.
   * @post  If !has_edge(_e_.node1(), _e_.node2()), then num_edges() remains the same.
   * @post  All pre-existing EdgeIterators are invalidated by this operation.
   * @post  All pre-existing IncidentIterators are invalidated by this operation.
   *
   * The worst-case complexity of this operation is
   * O(_e_.node1().degree() + _e_.node2().degree()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  /** Remove an edge from this graph.
   *
   * @param[in] e_it   An EdgeIterator corresponding to the edge to remove.
   *
   * @return    edge_begin(), valid for this graph after edge removal.
   *
   * @pre   e_it is a valid EdgeIterator for this graph.
   * @pre   e_it != edge_end().
   * @post  num_edges() is reduced by one.
   * @post  All pre-existing EdgeIterators are invalidated by this operation.
   * @post  All pre-existing IncidentIterators are invalidated by this operation.
   *
   * The worst-case complexity of this operation is
   * O((*e_it).node1().degree() + (*e_it).node2().degree()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }


private:
  // Internal type for nodes
  struct internal_node {
    Point point;            // position of this node in 3-d space
    node_value_type value;  // value of this node
    size_type index;        // index of this node in the graph
    internal_node(const Point& point, const node_value_type& value, size_type index) :
      point(point), value(value), index(index) {}
  };

  // Container for all nodes ever added to this graph (including those that may
  // have been removed); as a result, any node's uid_ is just an index into
  // this vector, which ensures fast access operations
  std::vector<internal_node> nodes_;

  // Container for uids of valid nodes in this graph, such that valid_uids_[i]
  // is the i-th node's uid, i.e., its index in nodes_, so that
  //    nodes_[valid_uids_[i]].index = i;
  std::vector<size_type> valid_uids_;

  // Number of edges in this graph
  size_type num_edges_;

  // Containers for uids of adjacent nodes to each node ever added to this
  // graph; the outer index is the uid of a node, while each value in the
  // corresponding sub-container is a uid for an adjacent node. Each of
  // adj_uids_max_ and adj_uids_min_ individually encode all edges in the graph
  // without repetition; having both of them separated allows us to:
  //  1) count up to the i-th edge without double-counting edges (see edge())
  //  2) return the degree of a node in O(1) time (see degree())

  // adj_uids_max_ encodes so-called "max edges" where "node1 < node2", i.e.,
  // uid < adj_uids_max_[uid][i] always holds (see add_edge())
  std::vector<std::vector<size_type>> adj_uids_max_;

  // adj_uids_min_ encodes so-called "min edges" where "node1 > node2", i.e.,
  // uid > adj_uids_min_[uid][i] always holds (see add_edge())
  std::vector<std::vector<size_type>> adj_uids_min_;

  // Container for values of edges in this graph; the outer container is indexed
  // by the uid of any node ever added to this graph, while the inner container
  // contains the values of incident "max edges" only to avoid repetition, in the
  // same order as adj_uids_max_[uid]
  std::vector<std::vector<edge_value_type>> edge_values_;
};

#endif // CME212_GRAPH_HPP
