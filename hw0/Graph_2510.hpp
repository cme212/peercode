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
   unsigned size_;
   unsigned num_edges_;
   std::vector<Point> positions_;
   // adjacency_matrix_[i][j] = 0 if there is no edge between i and j.
   // adjacency_matrix_[i][j] = idx if there is an edge.
   // idx is an index into edges_ map.
   // edges_ index starts with 1.
   // but edges are accessed using 0 index.
   std::vector<std::vector<unsigned>> adjacency_matrix_;
   std::map<int, std::pair<unsigned, unsigned>> edges_;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    size_ = 0;
    num_edges_ = 0;
    positions_ = std::vector<Point>();
    adjacency_matrix_ = std::vector<std::vector<unsigned>>();
    edges_ = std::map<int, std::pair<unsigned, unsigned>>();
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
  class Node {
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
    Node() {  }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->positions_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert (index_ < graph_->size_);
      return index_;
    }

    /** Return this node's graph reference. */
    Graph* graph() const {
      return graph_;
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
      bool graph_eq = (graph_ == n.graph());
      bool indx_eq = (index_ == n.index());
      return indx_eq && graph_eq;
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
      bool graph_eq = (graph_ == n.graph());
      return graph_eq && (index_ < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type index_;
    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *size_t
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    positions_.push_back(position);
    size_++;
    return Node(this, size_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    bool same_graph = (n.graph() == this);
    bool in_range = (n.index() < size_);
    return same_graph && in_range;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);        // Invalid node
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return this node's graph reference. */
    Graph* graph() const {
      return graph_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index_node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index_node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool same_graph = (e.graph() == e.graph_);
      bool same_edge = false;
      if (index_node1_ == e.node1().index() &&
          index_node2_ == e.node2().index()) {
        same_edge = true;
      }
      else if (index_node2_ == e.node1().index() &&
        index_node1_ == e.node2().index()) {
        same_edge = true;
      }
      return same_edge && same_graph;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      bool same_graph = (e.graph() == e.graph_);
      if (index_node1_ < e.node1().index() &&
          (index_node2_ < e.node2().index()))
        return true && same_graph;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type index_node1_;
    size_type index_node2_;
    size_type index_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type index_node1,
      size_type index_node2, size_type index)
        : graph_(const_cast<Graph*>(graph)) {
          index_node1_ = index_node1;
          index_node2_ = index_node2;
          index_ = index;
        }
  };

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
  Edge edge(unsigned i) const {
    // edges start indexing with 1 under the hood
    i++;
    return Edge(this,
      edges_.at(i).first,
      edges_.at(i).second,
      i);
  }

  std::pair<unsigned, unsigned> order_nodes(const Node& a, const Node& b) const
  {
    auto a_index = a.index();
    auto b_index = b.index();
    if (a.index() > b.index()) {
      a_index = b.index();
      b_index = a.index();
    }
    return std::make_pair(a_index, b_index);
  }
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  unsigned has_edge(const Node& a, const Node& b) const {
    auto ordered_pair = order_nodes(a, b);
    if (adjacency_matrix_.size() <= ordered_pair.first)
      return 0;
    if (adjacency_matrix_[ordered_pair.first].size() <= ordered_pair.second)
      return 0;
    return adjacency_matrix_[ordered_pair.first][ordered_pair.second];
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
    auto edge_check = has_edge(a, b);
    std :: cerr << (edge_check) << std :: endl ;
    if (edge_check > 0) {
      return edge(edge_check - 1);
    }
    else {
      auto ordered_pair = order_nodes(a, b);
      if (adjacency_matrix_.size() <= ordered_pair.first)
        adjacency_matrix_.resize(ordered_pair.first + 1,
          std::vector<unsigned>(ordered_pair.second + 1, 0));
      if (adjacency_matrix_[ordered_pair.first].size() <= ordered_pair.second)
          adjacency_matrix_[ordered_pair.first].resize(
            ordered_pair.second + 1, 0);
      num_edges_++;
      adjacency_matrix_[ordered_pair.first][ordered_pair.second] = num_edges_;
      edges_[num_edges_] = std::make_pair(a.index(), b.index());
      return Edge(this, a.index(), b.index(), num_edges_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    size_ = 0;
    num_edges_ = 0;
    positions_.clear();
    adjacency_matrix_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

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

  //
  // Incident Iterator
  //

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

  //
  // Edge Iterator
  //

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
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
