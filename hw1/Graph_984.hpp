#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>

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
  // Predeclare the internal structs for nodes and edges (as in proxy_example.cpp)
  struct internal_node;
  struct internal_edge;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Node value type. */
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
  Graph() : graph_nodes(), graph_edges(), next_node_id(0),
             next_edge_id(0), edge_nodes() {
  }

  /** Default destructor */
  ~Graph() = default; // Will generate default destructor

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
    }

    /** Find the internal_node element of this node. */
    internal_node& fetch() const {
      size_type this_index = graph_->node_id_index[id_];
      return graph_->graph_nodes[this_index];
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().location;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->node_id_index.at(id_);
    }

    /** Return node value. */
    node_value_type& value() {
      return fetch().value;
    }

    /** Return node value. */
    const node_value_type& value() const {
      return fetch().value;
    }

    /** Return degree of node. */
    size_type degree() const {
      return fetch().adj_nodes.size();
    }

    /** Return start incident iterator for this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this->graph_, 0, this->id_);
    }

    /** Return end incident iterator for this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, this->degree(), this->id_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.graph_ == this->graph_ && n.id_ == this->id_) {
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

     * A node is less than @a n if i) it comes from the same graph and has
     * a lower identification number, or ii) its graph < the other graph.
     */
    bool operator<(const Node& n) const {
      if (n.graph_ == this->graph_ && n.id_ < this->id_) {
        return true;
      }
      else if (n.graph_ < this->graph_) {
        return true;
      }
      return false;
    }

   private:
    // Private attributes: graph pointer and unique identifier
    graph_type* graph_;
    size_type id_;

    // Private constructor of valid Node objects
    Node(const graph_type* graph, size_type idx)
        : graph_(const_cast<graph_type*>(graph)), id_(idx) {
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return graph_nodes.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Create internal node
    internal_node new_element;
    new_element.location = position;
    new_element.id = next_node_id;
    new_element.value = value;

    // Add identifier-index combination to map
    next_node_id += 1;
    node_id_index[new_element.id] = this->size();

    // Return added node
    this->graph_nodes.push_back(new_element);
    return Node(this, new_element.id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_==this
         && this->node_id_index.find(n.id_)!=this->node_id_index.end()) {
        return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Use the vector of internal nodes to map index to identifier.
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, this->graph_nodes.at(i).id);
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->graph_edges[id_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->graph_edges[id_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     *
     * Edges are equal if they come from the same graph and have the same id.
     */
    bool operator==(const Edge& e) const {
      if (e.graph_==this->graph_ && e.id_==this->id_) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * An edge is less than another edge if i) both edges come from the
     * same graph and the edge has a smaller identifier, or ii) the edges
     * are from different graphs where (graph < other edge's graph).
     */
    bool operator<(const Edge& e) const {
      if (e.graph_==this->graph_ && e.id_ < this->id_) {
        return true;
      }
      else if (e.graph_ < this->graph_) {
        return true;
      }
      return false;
    }

   private:
    // Create edge attributes: graph pointer and edge identifier
    graph_type* graph_;
    size_type id_;

    // Constructor allowing the Graph class to construct valid edges
    Edge(const graph_type* graph, size_type idx)
       : graph_(const_cast<graph_type*>(graph)), id_(idx) {
    }

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return graph_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, this->graph_edges.at(i).id);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Order nodes by id (start_node has smallest id)
    Node start_node;
    Node end_node;
    if (a.id_ < b.id_) {
      start_node = a;
      end_node = b;
    }
    else {
      start_node = b;
      end_node = a;
    }

    // Create key for unordered map by concatenating identifiers in string
    std::string start_node_id = std::to_string(start_node.id_);
    std::string end_node_id = std::to_string(end_node.id_);
    std::string separator = ",";
    std::string map_key = start_node_id + separator + end_node_id;

    // Check if key in unordered map -> O(1)
    if (edge_nodes.find(map_key) != edge_nodes.end()) {
      return true;
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
    // Order nodes by id (start_node has smallest id)
    Node start_node;
    Node end_node;
    if (a.id_ < b.id_) {
      start_node = a;
      end_node = b;
    }
    else {
      start_node = b;
      end_node = a;
    }

    // Create key for unordered map by concatenating identifiers in string
    std::string start_node_id = std::to_string(start_node.id_);
    std::string end_node_id = std::to_string(end_node.id_);
    std::string separator = ",";
    std::string map_key = start_node_id + separator + end_node_id;

    // If key in unordered map (edge already exists), return existing edge
    // but do change order of node1() and node2() if necessary
    if (edge_nodes.find(map_key) != edge_nodes.end()) {
      size_type current_edge_id = edge_nodes[map_key];
      size_type current_edge_idx = edge_id_index[current_edge_id];
      size_type& current_node1_id  = graph_edges[current_edge_idx].node1.id_;
      size_type& current_node2_id = graph_edges[current_edge_idx].node2.id_;

      current_node1_id = a.id_;
      current_node2_id = b.id_;

      return Edge(this, edge_nodes[map_key]);
    }

    // Else, create new edge and add to graph
    internal_edge new_edge;
    new_edge.node1 = a;
    new_edge.node2 = b;
    new_edge.id = next_edge_id;

    // Add node b as adjacent node for node a and the other way around
    (a.fetch().adj_nodes).push_back(b.id_);
    (b.fetch().adj_nodes).push_back(a.id_);

    // Add edge and its nodes to unordered map
    edge_nodes[map_key] = new_edge.id;

    // Set index of identifier and set next identifier
    edge_id_index[new_edge.id] = num_edges();
    next_edge_id += 1;

    // Return added edge
    this->graph_edges.push_back(new_edge);
    return Edge(this, new_edge.id);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    graph_nodes.clear();
    graph_edges.clear();
  }

  //
  // Node Iterator
  //

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
    NodeIterator() {
    }

    /* Dereferencing operator. Returns current node iterator is at. */
    Node operator*() const {
      size_type id = this->graph_->graph_nodes[this->idx_].id;
      return Node(this->graph_, id);
    }

    /* Increment operator. */
    NodeIterator& operator++() {
      ++idx_;
      return (*this);
    }

    /* Equality operator. Checks whether the current node iterator is
       equal to a given node iterator. */
    bool operator==(const NodeIterator& node_iterator) const {
      if ((this->graph_ == node_iterator.graph_) && (this->idx_ == node_iterator.idx_)) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    graph_type* graph_;
    size_type idx_;

    /* Construct a valid NodeIterator */
    NodeIterator(const graph_type* graph, size_type idx) :
                     graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }
  };

  /** Set begin node iterator. */
  node_iterator node_begin() {
    return NodeIterator(this, 0);
  }

  /** Set end node iterator. */
  node_iterator node_end() {
    return NodeIterator(this, this->num_nodes());
  }

  //
  // Incident Iterator
  //

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
    IncidentIterator() {
    }

    /** Dereferencing operator. */
    Edge operator*() const {
      // Retrieve index and id of the nodes for this edge
      size_type this_node_idx = graph_->node_id_index[node_id_];
      size_type adj_node_id = graph_->graph_nodes[this_node_idx].adj_nodes[adj_idx_];
      size_type adj_node_idx = graph_->node_id_index[adj_node_id];

      // Set up key of edge in unordered map (string concatenation)
      std::string key_node1 = std::to_string(this->node_id_);
      std::string separator = ",";
      std::string key_node2 = std::to_string(adj_node_id);
      std::string key = key_node1 + separator + key_node2;

      // Make sure current node is node1 in iterator
      size_type edge_id = this->graph_->edge_nodes[key];
      size_type edge_idx = this->graph_->edge_id_index[edge_idx];
      graph_->graph_edges[edge_idx].node1 = graph_->node(this_node_idx);
      graph_->graph_edges[edge_idx].node2 = graph_->node(adj_node_idx);

      // Return edge
      return graph_->edge(edge_idx);
    }

    /** Increment operator. */
    IncidentIterator& operator++() {
      this->adj_idx_ = this->adj_idx_ + 1;
      return (*this);
    }

    /** Equality operator. */
    bool operator==(const IncidentIterator& inc_iterator) const {
      if ((this->graph_ == inc_iterator.graph_) && (this->adj_idx_ == inc_iterator.adj_idx_)
          && (this->node_id_ == inc_iterator.node_id_)) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    // Set class attributes
    graph_type* graph_;
    size_type adj_idx_;
    size_type node_id_;

    /* Construct a valid incident iterator. */
    IncidentIterator(const graph_type* graph, size_type adj_idx, size_type node_id) :
              graph_(const_cast<graph_type*>(graph)), adj_idx_(adj_idx), node_id_(node_id) {
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
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereferencing operator. */
    Edge operator*() {
      return this->graph_->edge(this->idx_);
    }

    /** Increment operator. */
    EdgeIterator& operator++() {
      this->idx_ = this->idx_ + 1;
      return (*this);
    }

    /** Equality operator. */
    bool operator==(const EdgeIterator& edge_iterator) const {
      if ((this->graph_ == edge_iterator.graph_) && (this->idx_ == edge_iterator.idx_)) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    // Set class attributes
    graph_type* graph_;
    size_type idx_;

    /* Construct a valid edge iterator. */
    EdgeIterator(const graph_type* graph, size_type idx) :
                 graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }
  };

  /** Set begin of edge iterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Set end of edge iterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  // Create struct for internal nodes of graph (as in proxy example)
  struct internal_node {
    size_type id;
    Point location; // Node proxy does not contain Point, this struct does
    node_value_type value;
    std::vector<size_type> adj_nodes;
  };

  // Create struct for internal edges of graph (as in proxy example)
  struct internal_edge {
    size_type id;
    Node node1;
    Node node2;
  };

  // Create vectors of nodes and edges in this graph
  std::vector<internal_node> graph_nodes;
  std::vector<internal_edge> graph_edges;

  // Declare identifiers of the next to-be-added node and edge
  size_type next_node_id;
  size_type next_edge_id;

  // Create unordered maps to match node and edge id's to indexes
  std::unordered_map<size_type, size_type> node_id_index;
  std::unordered_map<size_type, size_type> edge_id_index;

  // Create unordered map to link an edge's nodes to its id (fast lookup)
  std::unordered_map<std::string, size_type> edge_nodes;
};

#endif // CME212_GRAPH_HPP
