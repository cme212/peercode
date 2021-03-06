#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V> class Graph {
// Declaring size_type before private attributes
 public:
  using node_value_type = V;
  using size_type = unsigned;
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

  using vector_iterator_ = std::vector<size_type >::iterator;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Stores the nodes by IDs
  std::vector<std::pair<Point, node_value_type> > node_vec;
  //std::map<size_type, std::pair<size_type,size_type>> edge_map;
  std::vector<std::map<size_type, size_type> > adj_lookup;
  std::vector<std::vector<size_type> > adj_vec;
  std::vector<std::pair<size_type, size_type> > edge_ids;
  //std::map<std::pair<size_type, size_type>, size_type> edge_map;
  

 public:

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
 //using size_type = unsigned;

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
    
    Node() {
      // HW0: YOUR CODE HERE
      
    }

    Node(const Graph* graph, Graph::size_type id){
      this->graph = const_cast<Graph*>(graph);
      this->node_id = id;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      return graph->node_vec[node_id].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    Graph::size_type index() const {
      // HW0: YOUR CODE HERE
      return node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      return graph->node_vec[node_id].second;
    }

    const node_value_type& value() const{
      return graph->node_vec[node_id].second;
    }
    
    size_type degree() const{
      return graph->adj_vec[node_id].size();
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(graph, node_id, 0);
    }

    incident_iterator edge_end() const{
      size_type v_size = graph->adj_vec[node_id].size();
      return IncidentIterator(graph, node_id, v_size);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //assert(graph == n.graph);
      return (node_id == n.node_id) && (graph == n.graph);
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
      if (graph < n.graph) return true;
      if ((node_id < n.node_id) && (graph==n.graph)) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to graph node is in, and node id
    Graph* graph;
    size_type node_id;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_vec.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return node_vec.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    // HW0: YOUR CODE HERE
    size_type new_node_ind = num_nodes();
    node_vec.push_back(std::make_pair(position, val));
    // Adding new positioning for edge lookup
    std::map<size_type, size_type> newmap;

    // Adding entreis to adjacency lookup and map
    adj_lookup.push_back(newmap);
    adj_vec.push_back(std::vector<size_type>());

    return Node(this, new_node_ind);
    }      

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const{
    // HW0: YOUR CODE HERE
    return (n.graph == this) && (n.node_id < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert (i <= num_nodes());
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    // Constructor for edges
    Edge(const Graph* graph, Graph::size_type id, 
           std::pair<Graph::size_type, Graph::size_type> node_pair){
      this->graph = const_cast<Graph*>(graph);
      this->edge_id = id;
      this->node_pair = node_pair;
    }

    /** Return a node of this Edge */
    Graph::Node node1() const {
      // HW0: YOUR CODE HERE
      return graph->node(node_pair.first);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Graph::Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(node_pair.second);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //assert(graph == e.graph);
      //HW0: YOUR CODE HERE
      
      return (graph == e.graph) && (node_pair == e.node_pair);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph < e.graph) return true;
      if (graph > e.graph) return false;
      if (node_pair.first < e.node_pair.first) {
        return true;
      } else if (node_pair.first > node_pair.first){
        return false;
      } else {
          if (node_pair.second < e.node_pair.second) {
          return true;
        } else {
          return false;
        }
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
    Graph* graph;
    Graph::size_type edge_id;
    std::pair<Graph::size_type, Graph::size_type> node_pair;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE

    return edge_ids.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < num_edges());
    
    return Edge(this, i, edge_ids[i]);
    return Edge(); // Invalid edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph == this);
    assert(b.graph == this);

    auto a_ind = a.index();
    auto b_ind = b.index();

    assert(a_ind != b_ind);
    assert(a_ind < num_nodes());
    assert(b_ind < num_nodes());

    if (adj_lookup[a_ind].count(b_ind) > 0) return true;
    //if (adj_lookup[b_ind].count(a_ind) > 0) return true;

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

    assert(a.graph == this);
    assert(b.graph == this);
    auto a_ind = a.index();
    auto b_ind = b.index();
    assert(a_ind != b_ind);
    assert(a_ind < num_nodes());
    assert(b_ind < num_nodes());

    auto node_pair = std::make_pair(a_ind, b_ind);
    Graph::size_type new_edge_id = num_edges();

    // If has edge (check both directions), return that edge

    if (has_edge(a,b)){
      return Edge(this, adj_lookup[a_ind][b_ind], node_pair);
    }

    adj_lookup[a_ind][b_ind] = new_edge_id;
    adj_lookup[b_ind][a_ind] = new_edge_id;

    adj_vec[a_ind].push_back(b_ind);
    adj_vec[b_ind].push_back(a_ind);

    // Add nodes to edge id lookup
    edge_ids.push_back(node_pair);
    return Edge(this, new_edge_id, node_pair);    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //HW0: Your code here
    node_vec.clear();
    adj_lookup.clear();
    adj_vec.clear();
    edge_ids.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    NodeIterator(const Graph* graph, size_type node_id){
      this->graph = const_cast<Graph*>(graph);
      this->node_id = node_id;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return graph->node(node_id);
    }
    
    NodeIterator& operator++(){
      node_id++;
      return *this;
    }
    
    bool operator==(const NodeIterator& ni) const{
      return (node_id == ni.node_id) && (graph == ni.graph);
    }

    // bool operator != (const NodeIterator& ni) const{
    //   return !(this == ni);
    // }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph;
    size_type node_id;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
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
  class IncidentIterator : private equality_comparable<IncidentIterator>  {
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

    IncidentIterator(const Graph* graph, size_type node_id, size_type vec_ind){
      this->graph = const_cast<Graph*>(graph);
      this->node1_id = node_id;
      this->vector_ind = vec_ind;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    Edge operator*() const {
      std::vector<size_type> v = graph->adj_vec[node1_id];
      size_type node2_id = v[vector_ind];
      auto edge_ind = graph->adj_lookup[node1_id][node2_id];
      return Edge(graph, edge_ind, std::make_pair(node1_id, node2_id));
    }

    IncidentIterator& operator++(){
      vector_ind++;
      return *this;
    }
    bool operator==(const IncidentIterator& oii) const{
      return ((graph == oii.graph) && (node1_id == oii.node1_id) && \
              (vector_ind == oii.vector_ind));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph;
    size_type node1_id;
    //vector_iterator_ node2_iter;
    size_type vector_ind;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator> {
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

    EdgeIterator(const Graph* graph, size_type edge_id){
      this->graph = const_cast<Graph*>(graph);
      this->edge_id = edge_id;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return graph->edge(edge_id);
    }
    
    EdgeIterator& operator++() {
      edge_id ++;
      return *this;
    }
    
    bool operator==(const EdgeIterator& oi) const {
      return (graph == oi.graph) && (edge_id == oi.edge_id);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph;
    size_type edge_id;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  private:

  std::pair<size_type, size_type> create_node_pair(size_type n1, size_type n2) const{
    std::pair<size_type, size_type> node_pair;
    if (n1 < n2){
      node_pair.first = n1;
      node_pair.second = n2;
    } else {
      node_pair.first = n2;
      node_pair.second = n1;
    }
    return node_pair;
  }

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth


};

#endif //CME212_GRAPH_HPP
