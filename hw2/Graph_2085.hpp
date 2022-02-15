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
template <typename V, typename E>
class Graph {
 private:

  // Predeclare the internal structs
  struct internal_node;
  struct internal_edge;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  // HW1 3.2
  using node_value_type = V;

  // HW2 4.1
  using edge_value_type = E;

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
  Graph() : node_id_counter_(0), edge_id_counter_(0) {
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
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodeid_umap[fetch().nodeid];
    }

    /** Return modifiable node value. */
    node_value_type& value() {
      return fetch().value;
    }

    /** Return const node value. */
    const node_value_type& value() const {
      return fetch().value;
    }

    /** Return node's degree. */
    size_type degree() const {
      return fetch().degree;
    }

    /** Return an iterator to edges incident to the node. */
    incident_iterator edge_begin() {
      return IncidentIterator(graph_,*this,0);
    }

    /** One past the last incident edge. */
    incident_iterator edge_end() {
      return IncidentIterator(graph_,*this,degree());
    }

    /** Return modifiable node position. */
    Point& position() {
      return fetch().position;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
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
      if (graph_ == n.graph_ and
	  graph_->nodes_vec[this->index()].nodeid <
	  graph_->nodes_vec[n.index()].nodeid) {
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
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Node(const graph_type* graph, size_type nodeid)
        : graph_(const_cast<graph_type*>(graph)), id_(nodeid) {
    }
    /** Fetches the internal node. 
     *  Complexity: O(1).
     */
    internal_node& fetch() const {
      return graph_->nodes_vec[graph_->nodeid_umap.at(id_)];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
    
    // add internal node to nodes_vec
    size_type new_nodeid = node_id_counter_;
    node_id_counter_++;
    internal_node new_internal_node;
    new_internal_node.position = position;
    new_internal_node.value = val;
    new_internal_node.degree = 0;
    new_internal_node.nodeid = new_nodeid;
    nodes_vec.push_back(new_internal_node);

    // add node to incident umap, should add default empty vector as value
    incident_umap[new_nodeid];

    // add node to nodeid umap
    nodeid_umap[new_nodeid] = nodes_vec.size() - 1;
    
    return Node(this, new_nodeid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ((n.graph_ == this) and (n.index() < num_nodes()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, nodes_vec[i].nodeid);
  }

  /** Removes node and all edges incident to it.
   *  @pre 0 <= @a n1.index() < num_nodes()
   *  @post if graph has node: new num_nodes() = old num_nodes() - 1
   *        else: new num_nodes() = old num_nodes()
   *  @post g.node(i).index() == i for all i with 0 ≤ i < g.num_nodes()
   *  @ return 1 if has_node(n1) and 0 if !has_node(n1)
   *
   *  Complexity: O(deg(n).
   */
  size_type remove_node(const Node& n1) {
    if (!has_node(n1)) {
      return 0;
    }
    Node& n = const_cast<Node&>(n1);
    auto ii = n.edge_begin();
    for (; ii != n.edge_end();) {
      remove_edge(*ii);
      ii = n.edge_begin();
    }
    // update nodeid_umap
    nodeid_umap[nodes_vec.back().nodeid] = n.index();
    // swap and pop
    std::iter_swap(nodes_vec.begin() + n.index(), nodes_vec.end() - 1);
    nodes_vec.pop_back();
    return 1;
  }

  /** Removes node associated with node iterator and returns an iterator
   *  pointing to another node.
   *  @pre if graph has num_nodes > 0, node iterator points to a valid node
   *  @post if graph has node: new num_nodes() = old num_nodes() - 1
   *        else: new num_nodes() = old num_nodes()
   *  @post g.node(i).index() == i for all i with 0 ≤ i < g.num_nodes()
   *  @post if new num_nodes == 0, node_begin() == node_end()
   *  @return valid note_iterator
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
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
      return graph_->node(fetch().node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(fetch().node2_id);
    }

    /** Return modifiable edge value. */
    edge_value_type & value() {
      return fetch().value;
    }
    const edge_value_type& value() const {
      return fetch().value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
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
      if (graph_ == e.graph_ and id_ < e.id_) {
	return true;
      }
      else if (e.graph_ != graph_) {
	return true;
      }
      return false;
    }

    /** Return distance between node1 and node2.
     *  @post Node 1 and 2 position unchanged.
     *  Complexity: O(1).
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return distance between node1 and node2 when the edge was initialized.
     *  @post rest_length() unchanged.
     *  Complexity: O(1).
     */
    double rest_length() const {
      return fetch().rest_length;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type id_;

    Edge(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), id_(index) {}
    /** Fetches the internal edge from edges_vec.
     *  Complexity: O(1).
     */
    internal_edge& fetch() const {
      return graph_->edges_vec.at(graph_->edgeid_umap[id_]);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this,edges_vec[i].edgeid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::string edge_str = get_edge_str(a,b);
    return edges_umap.count(edge_str);
  }
  /** Add an edge to the graph, or return the current edge if it already exists
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
    if (has_edge(a,b) == true) {
      return Edge(this,get_edgeid(a,b));
    }
    // add internal edge to edge_vec
    size_type new_edgeid = edge_id_counter_;
    edge_id_counter_++;
    internal_edge new_internal_edge;
    new_internal_edge.node1_id = a.index();
    new_internal_edge.node2_id = b.index();
    new_internal_edge.edgeid = new_edgeid;
    new_internal_edge.rest_length = norm(a.position()-b.position());
    edges_vec.push_back(new_internal_edge);
    //update degree of both nodes
    nodes_vec[a.index()].degree++;
    nodes_vec[b.index()].degree++;
    // add edge to edge map
    std::string edge_str;
    if (a < b) {
      edge_str = std::to_string(nodes_vec[a.index()].nodeid) + ',' +
	std::to_string(nodes_vec[b.index()].nodeid);
    }
    else {
      edge_str = std::to_string(nodes_vec[b.index()].nodeid) + ',' +
	std::to_string(nodes_vec[a.index()].nodeid);
    }
    edges_umap[edge_str] = new_edgeid;
    // add edge to incident map
    incident_umap[nodes_vec[a.index()].nodeid].push_back(new_edgeid);
    incident_umap[nodes_vec[b.index()].nodeid].push_back(new_edgeid);
    
    // add edge to edgeid map
    edgeid_umap[new_edgeid] = edges_vec.size()-1;
    return Edge(this,new_edgeid);
  }

  /** Removes edge from graph.
   *  @pre n1 and n2 are valid nodes in the graph.
   *  @return 1 if graph has edge and 0 otherwise.
   *  @post If old has_edge(@a n1, @a n2), num_edges()--.
   *       Else,                        new num_edges() == old num_edges().
   *  @post Can invalidate edge indexes.
   *  @post has_edge(@a n1, @a n2) == false.
   *  Complexity: O(deg(n)).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1,n2)) {
      return 0;
    }
    size_type edgeid = get_edgeid(n1,n2);
    size_type edgeidx = edgeid_umap[edgeid];

    auto vec_iter_begin = incident_umap[nodes_vec[n1.index()].nodeid].begin();
    auto vec_iter_end = incident_umap[nodes_vec[n1.index()].nodeid].end();
    std::vector<size_type>::iterator it;
    it = find(vec_iter_begin,vec_iter_end,edgeid);
    incident_umap[nodes_vec[n1.index()].nodeid].erase(it);

    vec_iter_begin = incident_umap[nodes_vec[n2.index()].nodeid].begin();
    vec_iter_end = incident_umap[nodes_vec[n2.index()].nodeid].end();
    it = find(vec_iter_begin,vec_iter_end,edgeid);
    incident_umap[nodes_vec[n2.index()].nodeid].erase(it);
    
    std::string edge_str = get_edge_str(n1,n2);
    edges_umap.erase(edge_str);
    
    nodes_vec[n1.index()].degree--;
    nodes_vec[n2.index()].degree--;

    edgeid_umap[edges_vec.back().edgeid] = edgeidx;
    std::iter_swap(edges_vec.begin() + edgeidx, edges_vec.end() - 1);
    edges_vec.pop_back();
    edgeid_umap.erase(edgeid);
    return 1;
  }

  /** Removes edge from graph.
   *  @pre @a e is a valid edge in the graph
   *  Given node1 and node2 associated with edge @a e:
   *  @return 1 if graph has edge and 0 otherwise.
   *  @post If old has_edge(@a n1, @a n2), num_edges()--.
   *       Else,                        new num_edges() == old num_edges().
   *  @post Can invalidate edge indexes.
   *  @post has_edge(@a n1, @a n2) == false.
   *  Complexity: O(deg(n)).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),e.node2());
  }

  /** Removes edge associated with edge iterator and returns an iterator 
   *  pointing to another node.
   *  @pre if graph has num_edges > 0, iterator points to a valid edge
   *  @post if graph has edge: new num_edges() = old num_edges() - 1
   *        else: new num_edges() = old num_edges()
   *  @post Can invalidate node indexes.
   *  @post if new num_edges == 0, edge_begin() == edge_end()
   *  @return valid edge_iterator
   *  Complexity: O(deg(n)).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0.
   * Complexity: O(num_nodes() + num_edges()).
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeid_umap.clear();
    edgeid_umap.clear();
    incident_umap.clear();
    edges_umap.clear();
    nodes_vec.clear();
    edges_vec.clear();
    return;
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

    
    Node operator*() const {
      return node_val_;
    }
    NodeIterator& operator++() {
      idx_++;
      if (idx_ != graph_->nodes_vec.size()) {
	node_val_ = graph_->node(graph_->nodes_vec[idx_].nodeid);
      }
      return (*this);
    }
    bool operator==(const NodeIterator& node_itr) const {
      return (node_itr.graph_ == graph_ and
	      node_itr.idx_ == idx_);
    }
    
   private:
    friend class Graph;
    NodeIterator(const graph_type* graph_ptr, size_type idx)
      : graph_(const_cast<graph_type*>(graph_ptr)),
	idx_(idx)
    {
      if (idx_ != graph_->nodes_vec.size()) {
	node_val_ = graph_->node(graph_->nodes_vec[idx_].nodeid);
      }
    }
    graph_type* graph_;
    size_type idx_;
    value_type node_val_;
  };

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
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

    Edge operator*() const {
      return edge_val_;
    }
    IncidentIterator& operator++() {
      edge_incident_id_++;
      if (edge_incident_id_ != n_.degree()) {
	edge_id_ = graph_->incident_umap.at(nodeid_).at(edge_incident_id_);
	edge_val_ = graph_->edge(graph_->edgeid_umap.at(edge_id_));
	if (graph_->edges_vec.at(graph_->edgeid_umap.at(edge_id_)).node1_id !=
	    nodeid_) {
        // swap node1_id and node2_id (there may be cleaner code for this)
	  graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node2_id =
	    graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node1_id;
	  graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node1_id =
	    nodeid_;
	}
      }
      return (*this);
    }
    bool operator==(const IncidentIterator& inc_iter) const {
      return (inc_iter.edge_incident_id_ == edge_incident_id_);
    }
    
   private:
    friend class Graph;
    IncidentIterator(const graph_type* graph_ptr,
		     Node n,
		     size_type edge_incident_id)
      : graph_(const_cast<graph_type*>(graph_ptr)),
        n_(n),
	nodeid_(graph_ptr->nodes_vec[n.index()].nodeid),
        edge_incident_id_(edge_incident_id)
    {
      if (edge_incident_id_ != n_.degree()) {
	edge_id_ = graph_->incident_umap.at(nodeid_).at(edge_incident_id_);
	edge_val_ = graph_->edge(graph_->edgeid_umap[edge_id_]);
        if (graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node1_id !=
	    nodeid_) {
	  // swap node1_id and node2_id (there may be cleaner code for this)
	  graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node2_id =
	   graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node1_id;
	  graph_->edges_vec.at(graph_->edgeid_umap[edge_id_]).node1_id =
	    nodeid_;
	}
      }
    }
    graph_type* graph_;
    Node n_;
    size_type nodeid_;
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

    Edge operator*() const {
      return edge_val_;
    }
    EdgeIterator& operator++() {
      idx_++;
      if (idx_ != graph_->edges_vec.size()) {
	edge_val_ = graph_->edge(graph_->edges_vec[idx_].edgeid);
      }
      return (*this);
    }
    bool operator==(const EdgeIterator& edge_itr) const {
      return (edge_itr.graph_ == graph_ and
	      edge_itr.idx_ == idx_);
    }

   private:
    friend class Graph;
    EdgeIterator(const graph_type* graph_ptr, size_type idx)
      : graph_(const_cast<graph_type*>(graph_ptr)), idx_(idx) {
      if (idx_ != graph_->edges_vec.size()) {
        edge_val_ = graph_->edge(graph_->edges_vec[idx_].edgeid);
      }
    }
    graph_type* graph_;
    size_type idx_;
    value_type edge_val_;
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_vec.size());
  }

  
 private:

  size_type node_id_counter_;
  size_type edge_id_counter_;
  
  std::vector<internal_node> nodes_vec;
  std::vector<internal_edge> edges_vec;
  
  struct internal_node {
    Point position;
    node_value_type value;
    size_type degree;
    size_type nodeid;
  };
  struct internal_edge {
    size_type node1_id;
    size_type node2_id;
    size_type edgeid;
    edge_value_type value;
    double rest_length;
  };

  // nodeid map
  // nodeid -> index
  std::unordered_map<size_type,size_type> nodeid_umap;

  //edgeid map
  // edgeid -> index
  std::unordered_map<size_type,size_type> edgeid_umap;
  
  // node map
  // keys are node ids
  // values are a vector of edge ids incident to that node
  std::unordered_map<size_type,std::vector<size_type>> incident_umap;

  // map with pair of nodes as key and edge index as value
  std::unordered_map<std::string,size_type> edges_umap;

  /** Return edge id given two node ids.
   * @pre @a n1 and @a n2 are valid nodes in the graph.
   * @pre Graph has edge (@a n1, @a n2).
   * @return an edge id.
   * @post edges_umap unchanged
   *
   * Complexity: O(1), uses unordered_map, get_edge_str also O(1)
   */
  size_type get_edgeid(node_type n1, node_type n2) {
    std::string edge_str = get_edge_str(n1,n2);
    return edges_umap.at(edge_str);
  }

  /** Return a string to use as a key in edges_umap
   *  Concatenates the id's of the nodes with a comma.
   * @pre @a n1 and @a n2 are valid nodes in the graph.
   * @pre Graph has edge (@a n1, @a n2).
   * @post edge_str is a key present in edges_umap.
   * Complexity: O(1).
   */
  std::string get_edge_str(node_type n1, node_type n2) const {
    std::string edge_str;
    if (n1 < n2) {
      edge_str = std::to_string(nodes_vec[n1.index()].nodeid) + ',' +
        std::to_string(nodes_vec[n2.index()].nodeid);
    }
    else {
      edge_str = std::to_string(nodes_vec[n2.index()].nodeid) + ',' +
        std::to_string(nodes_vec[n1.index()].nodeid);
    }
    return edge_str;
  }
};

#endif // CME212_GRAPH_HPP
