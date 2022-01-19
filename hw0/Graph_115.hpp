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
class Graph {
private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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

  // ==========================================================================
  // CONSTRUCTORS AND DESTRUCTOR
  // ==========================================================================

  /** Construct an empty graph. */
  // HW0: YOUR CODE HERE
  Graph() : nodes_(), valid_uids_(), num_edges_(0), adj_uids_() {};

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
  class Node {
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
    // HW0: YOUR CODE HERE
    Node() : graph_(nullptr), uid_(0) {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(this->graph_ != nullptr);  // this node must belong to a graph
      return graph_->nodes_[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(this->graph_ != nullptr);  // this node must belong to a graph
      return graph_->nodes_[uid_].index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
      return std::tie(graph_, uid_) < std::tie(n.graph_, n.uid_);
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to this node's graph (can be nullptr if this node is invalid)
    graph_type* graph_;

    // Unique identifier used to access this node's data in the graph
    size_type uid_;

    /** Construct a valid Node.
     *  @param[in] graph  The graph containing this node.
     *  @param[in] uid    The unique identifier for this node in the graph.
     *  Complexity: O(1).
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
    // HW0: YOUR CODE HERE
    return valid_uids_.size();
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    nodes_.emplace_back(position, valid_uids_.size());
    valid_uids_.push_back(nodes_.size() - 1);
    adj_uids_.resize(nodes_.size());
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this                                // correct graph pointer
           && n.index() < n.graph_->size()                 // node index is in-bounds
           && n.uid_ == n.graph_->valid_uids_[n.index()];  // node uid matches the internal uid at the same index
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
  class Edge {
  public:
    /** Construct an invalid Edge. */
    // HW0: YOUR CODE HERE
    Edge() : graph_(nullptr), uid1_(0), uid2_(0) {}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return graph_ == e.graph_ && ((uid1_ == e.uid1_ && uid2_ == e.uid2_)
                                    || (uid1_ == e.uid2_ && uid2_ == e.uid1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      // We need to be careful of the case when comparing two edges with the
      // same nodes in opposite order; to preserve trichotomy, we first order
      // the nodes for an edge according to their uids before comparing
      size_type uid_min = std::min(uid1_, uid2_);
      size_type uid_max = std::max(uid1_, uid2_);
      size_type e_uid_min = std::min(e.uid1_, e.uid2_);
      size_type e_uid_max = std::max(e.uid1_, e.uid2_);
      return std::tie(graph_, uid_min, uid_max) < std::tie(e.graph_, e_uid_min, e_uid_max);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

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
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());  // avoid out-of-bounds issues
    size_type uid1;
    size_type adj_index = i;
    for (uid1 = 0; uid1 < adj_uids_.size(); ++uid1) {
      // If adj_index is too large to index into the vector of nodes adjacent
      // to the current node, decrement adj_index by the number of adjacencies
      // and continue to the next node
      if (adj_index >= adj_uids_[uid1].size()) {
        adj_index -= adj_uids_[uid1].size();
      } else {
      // Otherwise, we have found the node and the index of its adjacent node
      // that together form the i-th edge
        break;
      }
    }
    return Edge(this, uid1, adj_uids_[uid1][adj_index]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (a == b) {
      return false;  // no edges from a node to itself
    }
    size_type uid_min = std::min(a.uid_, b.uid_);
    size_type uid_max = std::max(a.uid_, b.uid_);
    for (size_type i = 0; i < adj_uids_[uid_min].size(); ++i) {
      if (adj_uids_[uid_min][i] == uid_max) return true;
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
    assert(!(a == b));  // no edges can be added from a node to itself
    if (!has_edge(a, b)) {
      // Record the larger uid of the two nodes forming an edge, at the index
      // of adj_uids_ matching the smaller uid
      size_type uid_min = std::min(a.uid_, b.uid_);
      size_type uid_max = std::max(a.uid_, b.uid_);
      adj_uids_[uid_min].push_back(uid_max);
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
    // HW0: YOUR CODE HERE
    nodes_.clear();
    valid_uids_.clear();
    num_edges_ = 0;
    adj_uids_.clear();
  }

  // ==========================================================================
  // NODE ITERATOR
  // ==========================================================================

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  // ==========================================================================
  // INCIDENT ITERATOR
  // ==========================================================================

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
  };

  // ==========================================================================
  // EDGE ITERATOR
  // ==========================================================================

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

  // Internal type for nodes
  struct internal_node {
    Point point;      // position of this node in 3-d space
    size_type index;  // index of this node in the graph
    internal_node(const Point& point, size_type index) :
      point(point), index(index) {}
  };

  // Container for all nodes ever added to this graph (including those that may
  // have been removed); as a result, any node's uid_ is just an index into
  // this vector, which ensures fast access operations
  std::vector<internal_node> nodes_;

  // Container for uids of valid nodes in this graph, such that valid_uids_[i]
  // is the i-th node's uid, i.e., its index in nodes_
  std::vector<size_type> valid_uids_;

  // Number of edges in this graph
  size_type num_edges_;

  // Container for uids of adjacent nodes to each node ever added to this
  // graph; the outer index is the uid of a node, while each value in the
  // corresponding sub-container is a uid for an adjacent node; this container
  // individually encodes all edges in the graph without repetition, and
  // uid < adj_uids_[uid][i] always holds (see add_edge)
  std::vector<std::vector<size_type>> adj_uids_;
};

#endif // CME212_GRAPH_HPP
