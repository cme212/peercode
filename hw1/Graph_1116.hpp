#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp;
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <unordered_map>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Predeclare the internal structs
  struct internal_node;
  struct internal_edge;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  // HW1 3.2
  using node_value_type = V;

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
    Node() {
      // HW0: YOUR CODE HERE
    }
    
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_vec.at(id_).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    node_value_type& value() {
      return graph_->nodes_vec.at(id_).value;
    }
    const node_value_type& value() const {
      return graph_->nodes_vec.at(id_).value;
    }
    size_type degree() const {
      return graph_->nodes_vec.at(id_).degree;
    }
    incident_iterator edge_begin() {
      return IncidentIterator(graph_,id_,0);
    }
    incident_iterator edge_end() {
      return IncidentIterator(graph_,id_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_ == n.graph_ and n.id_ == id_;
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
      if (graph_ == n.graph_ and id_ < n.id_) {
	return true;
      }
      else if (n.graph_ != graph_) {
	return true;
      }
      else {
	return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Node(const graph_type* graph, size_type nodeid)
        : graph_(const_cast<graph_type*>(graph)), id_(nodeid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_vec.size();
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
  Node add_node(const Point& position,
		const node_value_type& val = node_value_type ()) {
    // HW0: YOUR CODE HERE

    // add internal node to nodes_vec
    // THIS MIGHT NOT WORK IF NODES CAN BE REMOVED:
    // may set new_nodeid to a value already used and in the graph
    size_type new_nodeid = num_nodes();
    internal_node new_internal_node;
    new_internal_node.position = position;
    new_internal_node.value = val;
    new_internal_node.degree = 0;
    nodes_vec.push_back(new_internal_node);

    // add node to node map, should add default empty vector as value
    nodes_umap[new_nodeid];
    return Node(this, new_nodeid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // map.count returns 1 if in map, else 0
    return (n.graph_ == this) and (nodes_umap.count(n.id_));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0);
    assert(i < num_nodes());
    return Node(this, i);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(graph_->edges_vec[id_].node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(graph_->edges_vec[id_].node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ == e.graph_ and e.id_ == id_) {
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
      //HW0: YOUR CODE HERE
      if (graph_ == e.graph_ and e.id_ < id_) {
	return true;
      }
      else if (e.graph_ != graph_) {
	return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type id_;

    Edge(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), id_(index) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this,this->edges_vec.at(i).edgeid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(!(a==b));
    std::string edge_str = std::to_string(a.index()) + ',' + std::to_string(b.index());
    if (!(a < b)) {
	edge_str = std::to_string(b.index()) + ',' + std::to_string(a.index());
      }
    return edges_umap.count(edge_str);
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
    assert(!(a==b));
    if (has_edge(a,b) == true) {
      if (a < b) {
	return Edge(this,get_edgeid(a.index(),b.index()));
      }
      return Edge(this,get_edgeid(b.index(),a.index()));
    }
    // add internal edge to edge_vec
    size_type new_edgeid = num_edges();
    internal_edge new_internal_edge;
    new_internal_edge.node1_id = a.index();
    new_internal_edge.node2_id = b.index();
    new_internal_edge.edgeid = new_edgeid;
    edges_vec.push_back(new_internal_edge);
    //update degree of both nodes
    nodes_vec[a.index()].degree++;
    nodes_vec[b.index()].degree++;
    // add edge to edge map
    std::string edge_str;
    if (a < b) {
      edge_str = std::to_string(a.index()) + ',' + std::to_string(b.index());
    }
    else {
      edge_str = std::to_string(b.index()) + ',' + std::to_string(a.index());
    }
    edges_umap[edge_str] = new_edgeid;
    // add edge to node map
    nodes_umap[a.index()].push_back(new_edgeid);
    nodes_umap[b.index()].push_back(new_edgeid);
    return Edge(this,new_edgeid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    delete this;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    Node operator*() const {
      return node_val_;
    }
    NodeIterator& operator++() {
      //assert(node_ptr_ != node_end());
      node_id_++;
      if (node_id_ != graph_->nodes_vec.size()) {
	node_val_ = graph_->node(node_id_);
      }
      return (*this);
    }
    bool operator==(const NodeIterator& node_itr) const {
      return (node_itr.graph_ == graph_ and
	      node_itr.node_id_ == node_id_);
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(const graph_type* graph_ptr, size_type node_id)
      : graph_(const_cast<graph_type*>(graph_ptr)),
	node_id_(node_id)
    {
      if (node_id_ != graph_->nodes_vec.size()) {
	node_val_ = graph_->node(node_id_);
      }
    }
    graph_type* graph_;
    size_type node_id_;
    value_type node_val_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator(this, nodes_vec.size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
      return edge_val_;
    }
    IncidentIterator& operator++() {
      edge_incident_id_++;
      if (edge_incident_id_ != graph_->nodes_vec.at(node_id_).degree) {
	edge_id_ = graph_->nodes_umap.at(node_id_).at(edge_incident_id_);
	edge_val_ = graph_->edge(edge_id_);
      }
      return (*this);
    }
    bool operator==(const IncidentIterator& inc_iter) const {
      return inc_iter.edge_incident_id_ == edge_incident_id_;
    }
    
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const graph_type* graph_ptr,
		     size_type node_id,
		     size_type edge_incident_id)
      : graph_(const_cast<graph_type*>(graph_ptr)),
        node_id_(node_id),
        edge_incident_id_(edge_incident_id)
    {
      if (edge_incident_id_ != graph_->nodes_vec.at(node_id_).degree) {
	edge_id_ = graph_ptr->nodes_umap.at(node_id).at(edge_incident_id);
	edge_val_ = graph_->edge(edge_id_);
	if (graph_->edges_vec.at(edge_id_).node1_id != node_id_) {
        // swap node1_id and node2_id (there may be cleaner code for this)     
        graph_->edges_vec.at(edge_id_).node2_id =
          graph_->edges_vec.at(edge_id_).node1_id;
        graph_->edges_vec.at(edge_id_).node1_id = node_id_;
      }
      }
    }
    graph_type* graph_;
    size_type node_id_;
    size_type edge_incident_id_;
    size_type edge_id_;
    value_type edge_val_;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    Edge operator*() const {
      return edge_val_;
    }
    EdgeIterator& operator++() {
      edge_id_++;
      if (edge_id_ != graph_->edges_vec.size()) {
	edge_val_ = graph_->edge(edge_id_);
      }
      return (*this);
    }
    bool operator==(const EdgeIterator& edge_itr) const {
      return (edge_itr.graph_ == graph_ and
	      edge_itr.edge_id_ == edge_id_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const graph_type* graph_ptr, size_type edge_id)
      : graph_(const_cast<graph_type*>(graph_ptr)),
	edge_id_(edge_id)
    {
      if (edge_id_ != graph_->edges_vec.size()) {
        edge_val_ = graph_->edge(edge_id_);
      }
    }
    graph_type* graph_;
    size_type edge_id_;
    value_type edge_val_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_vec.size());
  }

  
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<internal_node> nodes_vec;
  std::vector<internal_edge> edges_vec;
  
  struct internal_node {
    Point position;
    node_value_type value;
    size_type degree;
  };
  struct internal_edge {
    size_type node1_id;
    size_type node2_id;
    size_type edgeid;
  };

  // node map
  // keys are node ids
  // values are a vector of edge ids incident to that node
  std::unordered_map<size_type,std::vector<size_type>> nodes_umap;

  // map with pair of nodes as key and edge index as value
  std::unordered_map<std::string,size_type> edges_umap;

  /** Return edge id given two node ids.
   * @pre @a n1_id and @a n2_id are in the graph and so is the edge (n1,n2).
   * @pre n1 < n2
   * @return an edge id.
   * @post edges_umap unchanged
   *
   * Complexity: O(1), uses unordered_map
   */
  size_type get_edgeid(size_type n1_id, size_type n2_id) {
    std::string edge_str = std::to_string(n1_id) + ',' + std::to_string(n2_id);
    return edges_umap.at(edge_str);
  }
};

#endif // CME212_GRAPH_HPP
