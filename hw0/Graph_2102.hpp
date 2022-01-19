#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

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

  // Declare struct for node's internal data
  struct internal_node;

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
    node_counter = 0;
    edge_counter = 0;

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
    Node() {
      is_valid_ = false;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return fetch().n_id;
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
      if ((this->n_id_ == n.index()) && (this->parent_graph_ == n.parent_graph_ )){ 
        return true;
      } else {
        return false;
      }
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
      if (this->n_id_ < n.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer to graph that node is in
    Graph* parent_graph_;
    // Unique id for node
    size_type n_id_; 
    // Flag for valid id
    bool is_valid_;
    
    /** Private Constructor for Graph, creates valid node*/
    Node(const Graph* parent_graph, size_type n_id) 
      : parent_graph_(const_cast<Graph*>(parent_graph)), n_id_(n_id){
        is_valid_= true;
    }
    
    // Helper method to return current internal node object
    internal_node& fetch() const {
      // return the internal_node that has the same index as n_id_
      // check if n_id_ is a key:
      if (parent_graph_->nodes_.count(n_id_)){
        return parent_graph_->nodes_.at(n_id_);
      } else {
        assert(false);
      }
    }

    // Helper method to return the parent graph that node is in
    Graph* current_graph() const {
      return parent_graph_;
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
    size_type new_id = node_counter;
    node_counter++;
    internal_node add_internal;
    add_internal.pos = position;
    add_internal.n_id = new_id;
    nodes_.insert({new_id, add_internal});
    Node add = Node(this, new_id);
    return add;        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    bool check_id = nodes_.count(n.index());
    if (check_id) {
      bool check_graph = (this == n.current_graph());
      if (check_graph) {
        return true;
      } 
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
    if (nodes_.count(i)){
      return Node(this, i);
    } else {
      assert(false);
    }
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
    Edge() {
      is_valid_ = false;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return an_edge_.first;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return an_edge_.second;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((an_edge_.first== e.node1()) && (an_edge_.second == e.node2())) {
        return true;
      } else if ((an_edge_.second == e.node1()) && (an_edge_.first == e.node2())) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (*this == e) {
        return false;
      } else if ((an_edge_.first < e.node1()) && (an_edge_.second < e.node2())) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // An edge has a pair of nodes
    std::pair<Node, Node> an_edge_;
    // An edge has an id associated with it
    size_type e_id_;
    // Flag to see if Edge is valid
    bool is_valid_;
    
    // Private constructor (used by Graph) creates valid edge
    Edge(const Node& a, const Node& b, size_type e_id){
      an_edge_.first = a;
      an_edge_.second = b;
      e_id_ = e_id;
      is_valid_ = true;
    }

    // Helper function that returns pair of nodes data type
    std::pair<Node, Node> get_edge() const {
      return an_edge_;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::pair<Node, Node> pair1(a, b);
    std::pair<Node, Node> pair2(b, a);
    // check pairs_ to see if pair1 or pair2 are keys
    bool check1 = pairs_.count(pair1);
    bool check2 = pairs_.count(pair2);
    return (check1 || check2);
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
    if (has_edge(a, b)) {
      std::pair<Node, Node> pair1(a, b);
      std::pair<Node, Node> pair2(b, a);
      // check pairs_ to see if pair1 or pair2 are keys
      bool check1 = pairs_.count(pair1);
      bool check2 = pairs_.count(pair2);
      if (check1) {
        size_type id = pairs_.at(pair1);
        return edges_.at(id);
      } else if (check2) {
        size_type id = pairs_.at(pair2);
        return edges_.at(id);
      } else {
        std::cout <<" something wrong with add_edge!!!" << std::endl;
      }
    }
    size_type edge_id = edge_counter;
    edge_counter++;

    Edge new_edge = Edge(a, b, edge_id);
    edges_.insert({edge_id, new_edge});
    pairs_.insert({{a,b}, edge_id});
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_counter = 0;
    edge_counter = 0;
    nodes_.clear();
    edges_.clear();
    pairs_.clear();
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

  // Define internal node structure 
  struct internal_node {
    Point pos;      // Position of a node
    size_type n_id; // Unique id for node
  };

  // map of nodes
  std::map<size_type, internal_node> nodes_;
  size_type node_counter;

  // map of edges
  std::map<size_type, Edge> edges_;
  size_type edge_counter;
  // map of pairs
  std::map<std::pair<Node,Node>, size_type> pairs_;
};

#endif // CME212_GRAPH_HPP
