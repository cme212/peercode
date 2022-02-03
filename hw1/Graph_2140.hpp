#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// Add the template
template<typename V>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // Internal node representation
  struct node_rep;
  // Internal edge representation
  struct edge_rep;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_map_type = const std::unordered_map<unsigned,node_rep>*;

  using edge_map_type = const std::unordered_map<unsigned, edge_rep>*;

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  using node_value_type = V;

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
  // Store nodes in a map, edges in a map, node_edge_check for adding edge
  Graph() 
      : node_map{}, node_size(0), edge_size(0), edge_map{} {
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
  class Node: private totally_ordered<Node> {
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
    // Public constructor, creates invalid node
    Node() {
    }

    /** Return this node's position. Checks the node map at index. */
    const Point& position() const {
      return graph_->node_map.at(ind_).p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return ind_;
    }

    /** Return this node's value of type V, allows for modifying value. */
    node_value_type& value() {
      return graph_->node_map.at(ind_).val;
    }

    /** Return this node's value of type V, doesn't allow for modification. */
    const node_value_type& value() const {
      return graph_->node_map.at(ind_).val;
    }

    /** Returns this node's degree, the number of adjacent nodes. */
    size_type degree() const {
      return graph_->node_map.at(ind_).adj_nodes.size();
    }

    /** Returns beginning incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this,graph_->node_map.at(ind_).
      adj_nodes.begin(),graph_->node_map.at(ind_).adj_nodes.end());
    }

    /** Returns end incident iterator 
    (corresponding to the last incident node)
    */
    incident_iterator edge_end() const {
      return IncidentIterator(this,graph_->node_map.at(ind_).adj_nodes.end(),
      graph_->node_map.at(ind_).adj_nodes.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Checks whether the graph and index are equal.
      if ((n.graph_ == graph_) && 
      (n.ind_ == ind_)) {
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
    bool operator<(const Node& n) const {
      if ((ind_ < n.ind_) || (n.graph_ != graph_ && ind_ == n.ind_)) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the graph
    Graph* graph_;
    // Node index
    size_type ind_;
     /** Private Constructor */
    Node(const Graph* graph, size_type ind)
      : graph_(const_cast<Graph*>(graph)), ind_(ind) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position, the node value
   * node value is set to a default if not supplied
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
  const node_value_type& val = node_value_type()) {
    // Inserts a node to the node_map with empty adjacency set to start
    std::unordered_set<size_type> empty_adj_set;
    node_map.insert({node_size,node_rep(position,node_size,val,
    empty_adj_set)});
    // Adds 1 to the node_size
    node_size += 1;
    // Returns node
    return Node(this,node_size-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Searches in the map.
    if (node_map.find(n.ind_) == node_map.end()) {
      return false;
    }
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
      return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Checks for undirected nature of the graph.
      if (((e.node1() == this->node1()) && (e.node2() == this->node2())) ||
      ((e.node2() == this->node1()) && (e.node1() == this->node2()))) {
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
      // Tests if the two nodes are less than each other.
      if ((node1_ < e.node1_) && (node2_ < e.node2_)) {
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Define an edge as containing two nodes.
    node_type node1_;
    node_type node2_;

    /* Private constructor */
    Edge(const node_type node1, const node_type node2)
      : node1_(node1), node2_(node2) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Return edge with two nodes at appropriate indices.
    return Edge(node(edge_map.at(i).node_ind1),node(edge_map.at(i).node_ind2));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Considers undirected nature of the graph
    if ((node_map.at(a.ind_).adj_nodes.find(b.ind_) 
    == node_map.at(a.ind_).adj_nodes.end()) &&
    (node_map.at(b.ind_).adj_nodes.find(a.ind_) 
    == node_map.at(b.ind_).adj_nodes.end())) {
      return false;
    }
    return true;
  }

  /** Add an edge to the graph, return the current edge if it already exists.
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
    // Check if the edge exists.
    if (has_edge(a,b)) {
      return Edge(a,b);
    }
    // Add to the edge map
    edge_map.insert({edge_size,edge_rep(a.ind_,b.ind_)});
    // Add to the node adjacency sets
    node_map.at(a.ind_).adj_nodes.insert(b.ind_);
    node_map.at(b.ind_).adj_nodes.insert(a.ind_);
    edge_size += 1;
    return Edge(a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_map.clear();
    edge_map.clear();
    node_size = 0;
    edge_size = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  //private totally_ordered<Node
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;// Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Dereferencing operator. 
    Returns the node corresponding to the graph at a given index.
    */
    Node operator*() const {
      return Node(graph_ptr,graph_ptr->node_map.at(iter_idx).node_idx);
    }
    
    /** Advances the iterator. */
    NodeIterator& operator++() {
      if (iter_idx < graph_ptr->num_nodes()) {
        ++iter_idx;
      }
      return *this;
    }

    /** Iterators are equal if graphs and node indices are equal. */
    bool operator==(const NodeIterator& node_iter) const {
      return node_iter.graph_ptr == graph_ptr 
      && node_iter.iter_idx == iter_idx;
    }


   private:
    friend class Graph;

    const graph_type* graph_ptr;
    size_type iter_idx;

    // Private constructor
    NodeIterator(const graph_type* graph_ptr_, size_type iter_idx_) :
    graph_ptr(graph_ptr_), iter_idx{iter_idx_} {
    }
  };

  /** Begin node iterator, node at 0. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** End node iterator, one past the end. */
  node_iterator node_end() const {
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
    using reference         = Edge&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Dereferencing operator, returns edge attached
    to initial node and current adjacent node.
    */
    Edge operator*() const {
      return Edge(*node_ptr,Node((*node_ptr).graph_,*it));
    }

    /** Increment operator. */
    IncidentIterator& operator++() {
      if (it != end) {
        ++it;
      }
      return *this;
    }

    /** Equality operator. */
    bool operator==(const IncidentIterator& inc_iter) const {
      return node_ptr == inc_iter.node_ptr 
      && it == inc_iter.it && end == inc_iter.end;
    }

   private:
    friend class Graph;

    // stores a node pointer, iterator to beginning of adjacency set
    // iterator to end of adjacency set
    const Node* node_ptr;
    std::unordered_set<size_type>::iterator it;
    std::unordered_set<size_type>::iterator end;   

    IncidentIterator(const Node* node_ptr_, 
    std::unordered_set<size_type>::iterator it_,
    std::unordered_set<size_type>::iterator end_) :
    node_ptr(node_ptr_), it{it_}, end{end_} {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy



    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Returns the two nodes connected to this edge. */
    Edge operator*() const {
      return Edge(Node(graph_ptr,graph_ptr->edge_map.at(iter_idx).node_ind1),
      Node(graph_ptr,graph_ptr->edge_map.at(iter_idx).node_ind2));
    }

    /** ++ operator. */
    EdgeIterator& operator++() {
      if (iter_idx < graph_ptr->num_edges()) {
        ++iter_idx;
      }
      return *this;
    }

    /** Equal if they have the same graph and edge index. */
    bool operator==(const EdgeIterator& edge_iter) const {
      return edge_iter.graph_ptr == graph_ptr && edge_iter.iter_idx == iter_idx;
    }

   private:
    friend class Graph;

    // Edge iterator has graph pointer and iteration index.
    const graph_type* graph_ptr;
    size_type iter_idx;

    EdgeIterator(const graph_type* graph_ptr_, size_type iter_idx_) :
    graph_ptr(graph_ptr_), iter_idx{iter_idx_} {
    }
  };

  /** Edge iterator begins at first index. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /** Edge iterator end. */
  edge_iterator edge_end() const {
    return EdgeIterator(this,this->num_edges());
  }


 private:

  // structure for internal nodes, contains a point and node index
  // value and set of adjacent nodes
  struct node_rep {
    const Point p;
    size_type node_idx;
    node_value_type val;
    std::unordered_set<size_type> adj_nodes;
    node_rep(const Point curr_p, size_type curr_idx, node_value_type v, 
    std::unordered_set<size_type> adj_nodes_) : p(curr_p), 
    node_idx(curr_idx), val(v), adj_nodes(adj_nodes_) {}
  };

  // structure for internal edge, contains two node indices
  struct edge_rep {
    size_type node_ind1;
    size_type node_ind2;
    edge_rep(size_type i1, size_type i2) : node_ind1(i1), node_ind2(i2) {}
  };

  std::unordered_map<size_type,node_rep> node_map;
  size_type node_size;
  size_type edge_size;
  std::unordered_map<size_type, edge_rep> edge_map;

};

#endif // CME212_GRAPH_HPP
