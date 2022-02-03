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

  /** Type of user-specified value to use as Node attribute */
  using node_value_type = V;

  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), node_values_(), points_(), adj_(), edge_ptrs_(), edges_() {
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
    Node()
      // HW0: YOUR CODE HERE
      : graph_(nullptr), index_(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->points_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return size_type(index_);
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Set this node's user-specified value */
    void set_value(const size_type& val) {
      graph_->node_values_[index_] = val;
      return;
    }
    
    /** Return this node's user-specified value */
    node_value_type& value() {
      return graph_->node_values_[index_];
    }

    const node_value_type& value() const {
      return graph_->node_values_[index_];
    }


    /** Return the number of incident edges to this node */
    size_type degree() const {
      return graph_->adj_[index_].size();
    }

    /** Return begin and end incident iterators for this node */
    incident_iterator edge_begin() const {
      incident_iterator begin(graph_, index_, 0);
      return begin;
    }

    /** Return begin and end incident iterators for this node */
    incident_iterator edge_end() const {
      incident_iterator end(graph_, index_, graph_->adj_[index_].size());
      return end;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_ && index_ == n.index_);
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
      return (index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    graph_type* graph_;
    // This node's unique index
    size_type index_;

    /** Private Constructor */
    Node(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Add node to end of vector of nodes
    nodes_.push_back(node_type(this, nodes_.size()));
    // Add point to end of vector of points
    points_.push_back(position);
    // Add node value to vector of user-specified node values
    node_values_.push_back(node_value);
    // Add empty vector to adjacency vector of vectors
    adj_.push_back({});

    return nodes_.back();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n == nodes_[n.index()]);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Make sure index i is valid
    assert(i < num_nodes());
    assert(i >= 0);

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
    Edge() 
      // HW0: YOUR CODE HERE
      : graph_(nullptr), index_(0), node1_(0), node2_(0) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (index_ == e.index_ && graph_ == e.graph_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return index_ < e.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    graph_type* graph_;
    // This edge's unique index
    size_type index_;
    // Indices of the nodes that make up this edge
    size_type node1_;
    size_type node2_;
    /** Private Constructor */
    Edge(const graph_type* graph, size_type index, size_type node1, size_type node2)
        : graph_(const_cast<graph_type*>(graph)), index_(index), node1_(node1), node2_(node2) {
    }
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
    // Make sure index i is valid
    assert(i >= 0);
    assert(i < num_edges());

    // Loop through map values (edge objects) and check index of each one
    for (auto& nodes_edge: edges_) {
          if (nodes_edge.second.index_ == i) {
              return nodes_edge.second;
          }
      }

    // This statement will never be reached, but must provide some return
    // value for all flows to avoid warning
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Check that a and b are valid nodes of this graph
    assert(this == a.graph_ && this == b.graph_);

    return edges_.count(std::make_tuple(a.index(), b.index())) +
           edges_.count(std::make_tuple(b.index(), a.index()));
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
    // Check that a and b are valid nodes of this graph
    assert(this == a.graph_ && this == b.graph_);

    /** Check if graph already contains edge. has_edge is not used because
     * need to know specifically if it contains (a,b) or (b,a), but has_edge
     * does not distinguish.
     */
    if (edges_.count(std::make_tuple(a.index(), b.index()))) {
        return edges_.at(std::make_tuple(a.index(), b.index()));
    }
    else if (edges_.count(std::make_tuple(b.index(), a.index()))) {
        return edges_.at(std::make_tuple(b.index(), a.index()));
    }
    /** Add edge if graph does already contain it. Also, modify adjacency
     * vector at both node index locations.
     */
    else {
        edges_.insert({std::make_tuple(a.index(), b.index()),
                      Edge(this, edges_.size(), a.index(), b.index())});
        adj_[a.index()].push_back(&(edges_.at(std::make_tuple(a.index(), b.index()))));
        adj_[b.index()].push_back(&(edges_.at(std::make_tuple(a.index(), b.index()))));
        edge_ptrs_.push_back(&(edges_.at(std::make_tuple(a.index(), b.index()))));
        return edges_.at(std::make_tuple(a.index(), b.index()));
    }
  }

  /** Get maximum user-specified node value */
  V max_value() {
    return *(std::max_element(node_values_.begin(), node_values_.end()));
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    points_.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() 
      : graph_(nullptr), index_(0) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    Node operator*() const {
        assert(index_ < graph_->nodes_.size());
        return graph_->nodes_[index_];
    }
    
    NodeIterator& operator++() {
        index_++;
        assert(index_ <= graph_->nodes_.size());
        return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
        return graph_ == node_iter.graph_ && index_ == node_iter.index_;
    }

    bool operator!=(const NodeIterator& node_iter) const {
        return graph_ != node_iter.graph_ || index_ != node_iter.index_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Pointer back to the Graph container
    graph_type* graph_;
    // Current index
    size_type index_;
    /** Private Constructor */
    NodeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
    node_iterator begin(this, 0);
    return begin;
  }

  node_iterator node_end() const {
    node_iterator end(this, nodes_.size());
    return end;
  }

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
    IncidentIterator() : 
      graph_(nullptr), node_index_(0), adj_index_(0) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const {
      return *(graph_->adj_[node_index_][adj_index_]);
    }

    IncidentIterator& operator++() {
      adj_index_++;
      assert(adj_index_ <= graph_->adj_[node_index_].size());
      return *this;
    }

    bool operator==(const IncidentIterator& iter) const {
      return graph_ == iter.graph_ && 
             node_index_ == iter.node_index_ &&
             adj_index_ == iter.adj_index_;
    }

    bool operator!=(const IncidentIterator& iter) const {
      return graph_ != iter.graph_ ||
             node_index_ != iter.node_index_ ||
             adj_index_ != iter.adj_index_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Pointer back to the Graph container
    graph_type* graph_;
    // Current index
    size_type node_index_;
    size_type adj_index_;
    /** Private Constructor */
    IncidentIterator(const graph_type* graph, size_type node_index, size_type adj_index)
        : graph_(const_cast<graph_type*>(graph)),
          node_index_(node_index),
          adj_index_(adj_index) {
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
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() 
      : graph_(nullptr), index_(0) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      assert(index_ < graph_->edges_.size());
      return *(graph_->edge_ptrs_[index_]);
    }

    EdgeIterator& operator++() {
      index_++;
      assert(index_ <= graph_->edges_.size());
      return *this;
    }

    bool operator==(const EdgeIterator& edge_iter) const {
      return graph_ == edge_iter.graph_ && index_ == edge_iter.graph_;
    }

    bool operator!=(const EdgeIterator& edge_iter) const {
      return graph_ != edge_iter.graph_ || index_ != edge_iter.index_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer back to the Graph container
    graph_type* graph_;
    // Current index
    size_type index_;
    /** Private constructor */
    EdgeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin() const {
    edge_iterator begin(this, 0);
    return begin;
  }

  edge_iterator edge_end() const {
    edge_iterator end(this, edges_.size());
    return end;
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<node_type> nodes_;
  std::vector<node_value_type> node_values_;
  std::vector<Point> points_;
  std::vector<std::vector<edge_type*>> adj_;
  std::vector<edge_type*> edge_ptrs_;
  std::map<std::tuple<size_type, size_type>, edge_type> edges_;

};

#endif // CME212_GRAPH_HPP
