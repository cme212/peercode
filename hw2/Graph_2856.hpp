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

template <typename V, typename E>
class Graph {

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

  /** Synonym for the type of Node value. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for the type of Edge value. */
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


  private:

  // Vector of pointers to Point objects,  storing all the nodes of the graph
  std::vector<Point*> nodes;

  // Vector of values of the nodes (values[i] = nodes[i].value())
  std::vector<node_value_type> node_values;

  // Vector of user facing indices for nodes
  std::vector<size_type> node_idx;

  // Vector of mapping from idx to uid (internal id) for nodes
  std::vector<size_type> node_i2u;

  // Adjacency list (adjacency_list[i] is a list of all the adjacent node indexes of node i)
  std::vector<std::vector<size_type>> adjacency_list;

  // Map indices to edges (useful to accelerate edge(i))
  std::vector<Edge> edges;

  // Vector of values of the edges (values[i] = edges[i].value())
  std::vector<edge_value_type> edge_values;

  // Vector of user facing indices for edges
  std::vector<size_type> edge_idx;

  // Vector of mapping from idx to uid (internal id) for edges
  std::vector<size_type> edge_i2u;

  // Number of edges
  size_type num_edge;

  


  public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(){
    this->nodes = std::vector<Point*>();
    this->node_values = std::vector<node_value_type>();
    this->node_idx = std::vector<size_type>();
    this->node_i2u = std::vector<size_type>();
    this->adjacency_list = std::vector<std::vector<size_type>>();
    this->edges = std::vector<Edge>();
    this->edge_values = std::vector<edge_value_type>();
    this->edge_idx = std::vector<size_type>();
    this->edge_i2u = std::vector<size_type>();
    this->num_edge = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  /** Removes a node in the graph
   * @param[in]     n           The node to be removed
   * @return        1 if any node has been removed, 0 otherwise
   *
   * @pre g.node(i).index() == i for all 0 <= i < g.num_nodes()
   * @pre g.node(n.index()) == n  
   * 
   * @post g.node(i).index() == i for all 0 <= i < g.num_nodes()
   * 
   * Complexity: No more than O(num_nodes()), hopefully less
   */
  size_type remove_node (const Node& n){

    if (!this->has_node(n)) {return 0;}

    // We first remove all connected edges
    auto begin = n.edge_begin();
    while (begin != n.edge_end()){
      remove_edge(*(n.edge_begin()));
    }

    // Update the internal to user node vector
    size_type idx = node_idx[n.node_index];
    node_i2u[idx] = node_i2u.back();
    node_i2u.pop_back();
    node_idx[node_i2u[idx]] = idx;

    return 1;
  }

  /** Removes a node in the graph
  * @param[in]     n_it   The node_iterator corresponding to the node to be removed
  * @return        The node_iterator corresponding to the node we removed
  */
  node_iterator remove_node (node_iterator n_it){
    this->remove_node(*n_it);
    return n_it;
  }

  /** Removes an edge in the graph
   * @param[in]     a   The first extremity of the edge to be removed
   * @param[in]     b   The second extremity of the edge to be removed
   * 
   * @return        1 if any edge has been removed, 0 otherwise
   *
   * @post num_edges() should return the number of unique undirected edges
   * 
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge (const Node& a, const Node& b){

    if(! this->has_edge(a, b)){return 0;}

    size_type i_a = a.node_index;
    size_type i_b = b.node_index;

    // Find this edge's index
    size_type index;
    for(size_type i = 0; i < (this->edges).size(); i++){
      // The index is not taken into acount in the == operator
      if(this->edges[i] == Edge(this, i_a, i_b, 0)){
        index = i;
        break;
      }
      if(this->edges[i] == Edge(this, i_b, i_a, 0)){
        index = i;
        break;
      }
    }
    
    // Update the adjacency list
    for (auto it = (adjacency_list[i_a]).begin(); it != (adjacency_list[i_a]).end(); ++it)
    {
      if(*it == i_b){
        (adjacency_list[i_a]).erase(it);
        break;
      }
    }

    for (auto it = (adjacency_list[i_b]).begin(); it != (adjacency_list[i_b]).end(); ++it)
    {
      if(*it == i_a){
        (adjacency_list[i_b]).erase(it);
        break;
      }
    }
    
    // Update the edge internal to user list and num_edge
    size_type idx = edge_idx[index];
    edge_i2u[idx] = edge_i2u.back();
    edge_i2u.pop_back();
    edge_idx[edge_i2u[idx]] = idx;
    num_edge --;
    
    return 1;

  }

  /** Removes an edge in the graph
  * @param[in]     e   The edge to be removed
  * @return        1 if any edge has been removed, 0 otherwise
  */
  size_type remove_edge (const Edge& e){
    return(this->remove_edge(e.node1(), e.node2()));
  }

  /** Removes an edge in the graph
  * @param[in]     e_it   The edge_iterator corresponding to the edge to be removed
  * @return        The edge_iterator corresponding to the edge we removed
  */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge e = *e_it;
    remove_edge(e);
    return e_it;
  }

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

    // Constructor of a Node
    Node(const Graph* graph, size_type node_index){
      this->graph = const_cast<Graph*>(graph);
      this->node_index = node_index;
    }

    // HW2
    /** Returns a reference to the modifiable node's position */
    Point& position(){
      return *(((*(this->graph)).nodes)[this->node_index]);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return *(((*(this->graph)).nodes)[this->node_index]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return ((*(this->graph)).node_idx[this->node_index]);
    }

      /* Get access to the value of the current node.
     * @return the value of the current node object
     */
    node_value_type& value(){
      size_type i = this->node_index;
      return((*graph).node_values[i]);
    }

      /* Get access to the value of the current node.
     * @return the value of the current node object
     */
    const node_value_type& value() const{
      size_type i = this->node_index;
      return((*graph).node_values[i]);
    }

      /* Get access to the degree of the current node.
     * @return the number of neighbors of the current node object
     */
    size_type degree() const{
      size_type i = this->node_index;
      return(((*(this->graph)).adjacency_list[i]).size());
    }

      /* Get access to the first iterator when we iterate on the current node's neighbors.
     * @return the the first iterator when we iterate on the current node's neighbors
     */
    incident_iterator edge_begin() const{
      return(IncidentIterator(this->graph, this->node_index, 0));
    }

      /* Get access to the last iterator when we iterate on the current node's neighbors.
     * @return the the last iterator when we iterate on the current node's neighbors
     */
    incident_iterator edge_end() const{
      return(IncidentIterator(this->graph, this->node_index, this->degree()));
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return(this->graph == n.graph && this->node_index == n.node_index);
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
      if(this->graph < n.graph){return true;}
      return(this->graph == n.graph && this->node_index < n.node_index);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph;
    size_type node_index;


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return (this->node_i2u).size();
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

  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {

    // HW0: YOUR CODE HERE
    Point* pos = new Point;
    *pos = position;

    // We add the new node to the end of our list of nodes
    (this->nodes).push_back(pos);
    (this->node_i2u).push_back(nodes.size()-1);
    (this->node_idx).push_back(size()-1);

    // We add the new node's value to the list of values
    (this->node_values).push_back(v);

    // Node that we return
    Node res = Node(this, (this->nodes).size()-1);

    // Add an empty list of neighbours in the adjacency list
    this->adjacency_list.push_back(std::vector<size_type>(0));

    return res;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return(n.graph == this && n.index() < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < (this->node_i2u).size());
    return Node(this, node_i2u[i]);
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


  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Construct an Edge. */
    Edge(const Graph* graph, size_type node1_index, size_type node2_index, 
        size_type edge_id){
      this->graph = const_cast<Graph*>(graph);
      this->node1_index = node1_index;
      this->node2_index = node2_index;
      this->edge_index = edge_id;
    }

    /** Return the value of this Edge */
    edge_value_type& value(){
      return((*graph).edge_values[this->edge_index]);
    }

    /** Return the value of this Edge */
    const edge_value_type& value() const{
      return((*graph).edge_values[this->edge_index]);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph, this->node1_index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph, this->node2_index);
    }

    /** Return the length of an edge */
    double length() const {
      return(norm((this->node1()).position() - (this->node2()).position()));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return((this->graph == e.graph) && ((this->node1() == e.node1() && this->node2() == e.node2()) ||
      (this->node1() == e.node2() && this->node2() == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE

      // Defining the min and max nodes of the current edge and e
      Node min_node = Node(this->graph, std::min(this->node1().node_index, this->node2().node_index));
      Node max_node = Node(this->graph, std::max(this->node1().node_index, this->node2().node_index));
      Node min_node_e = Node(this->graph, std::min(e.node1().node_index, e.node2().node_index));
      Node max_node_e = Node(this->graph, std::max(e.node1().node_index, e.node2().node_index));

      if(this->graph < e.graph){return true;}

      // First we compare minimal nodes
      else if(this->graph == e.graph && min_node<min_node_e){return true;}

      // Then we compare maximal nodes
      else if(this->graph == e.graph && min_node==min_node_e && max_node<max_node_e){
        return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Indices of the two nodes
    size_type node1_index;
    size_type node2_index;

    // Index of the edge
    size_type edge_index;

    // Graph to which the edge belongs
    Graph* graph;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->num_edge;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges[this->edge_i2u[i]];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    if(!this->has_node(a) || !this->has_node(b)){return false;}

    for (size_type i = 0; i < this->adjacency_list[a.node_index].size(); i++)
    {
      if (this->adjacency_list[a.node_index][i] == b.node_index) {return true;}
    }

    for (size_type i = 0; i < this->adjacency_list[b.node_index].size(); i++)
    {
      if (this->adjacency_list[b.node_index][i] == a.node_index) {return true;}
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

    // See the comment at the declaration of the adjacency_list
    if(this->has_edge(a, b)){
      // Find this edge's index
      unsigned int index;
      for(unsigned int i = 0; i < (this->edges).size(); i++){
        // The index is not taken into acount in the == operator
        if(this->edges[i] == Edge(this, a.node_index, b.node_index, 0)){
          index = i;
        }
        if(this->edges[i] == Edge(this, b.node_index, a.node_index, 0)){
          index = i;
        }
      }
      return Edge(this, a.node_index, b.node_index, index);
    }

    else{
      assert(this->has_node(a) && this->has_node(b));
    
      // Update the adjacency list
      this->adjacency_list[a.node_index].push_back(b.node_index);
      this->adjacency_list[b.node_index].push_back(a.node_index);

      // Increment the number of edges
      this->num_edge ++;

      // Update the list of edges
      this->edges.push_back(Edge(this, a.node_index, b.node_index, edges.size()));
      
      // Update the list of values
      this->edge_values.push_back(edge_value_type());

      // Update edge_i2u
      this->edge_i2u.push_back(edges.size()-1);

      // Update edge_idx
      this->edge_idx.push_back(num_edge-1);



      return Edge(this, a.node_index, b.node_index, edges.size()-1);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Since nodes is a vector of pointers, we first deallocate
    // the memory of the points to avoid memory leak
    unsigned int n = this->nodes.size();
    for (unsigned int i = 0; i < n; i++){delete this->nodes[i];}

    // We clear the vectors
    this->nodes.clear();
    this->node_values.clear();
    this->node_i2u.clear();
    this->node_idx.clear();
    this->adjacency_list.clear();
    this->edges.clear();
    this->edge_values.clear();
    this->edge_i2u.clear();
    this->edge_idx.clear();
    this->num_edge = 0;
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

      /* Dereference operator
     * @return the node behind the current iterator
     */
    Node operator*() const{
      return(Node(graph, (*graph).node_i2u[this->node_index]));
    }

      /* Increments to the next node in the graph
     * @return the node iterator pointing at the next node in the graph
     */
    NodeIterator& operator++(){
      // Update the index
      this->node_index ++;

      return(*this); 
    }

      /* Defines equality between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is equal to the argument, false otherwise
     */
    bool operator==(const NodeIterator& it) const{
      return(this->graph == it.graph && this->node_index == it.node_index);
    }

     /* Defines difference between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is different from the argument, false otherwise
     */
    bool operator!=(const NodeIterator& it) const{
      return(!(this->graph == it.graph && this->node_index == it.node_index));
    }

   private:
    friend class Graph;

    // Graph in which we are
    Graph* graph;

    // Index of the node at which we are in this graph
    size_type node_index;


    NodeIterator(const Graph* g, size_type node_id){
      this->graph = const_cast<Graph*>(g);
      this->node_index = node_id;
    }
  };

   /* Return an iterator pointing at the first node of the graph
  * @return an iterator pointing at the first node of the graph
  */
  node_iterator node_begin() const{
    return(NodeIterator(this, 0));
  }

   /* Return an iterator pointing at the last node of the graph
  * @return an iterator pointing at the larst node of the graph
  */
  node_iterator node_end() const{
    return(NodeIterator(this, size()));
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
    IncidentIterator() {
    }

      /* Dereference operator
     * @return the edge behind the current iterator
     */
    Edge operator*() const{

      size_type i = this->node_index;

      // We look for the current adjacent node in the adjacency list
      size_type adj_node_index = (*this->graph).adjacency_list[i][this->edge_adj_id];

      // Find the edge index
      size_type index = 0;
      for(size_type k = 0; k < ((this->graph)->edges).size(); k++){
        // The index is not taken into acount in the == operator
        if(((this->graph)->edges)[k] == Edge(this->graph, i, adj_node_index, 0)){
          index = k;
          break;
        }
        if(((this->graph)->edges)[k] == Edge(this->graph, adj_node_index, i, 0)){
          index = k;
          break;
        }
      }

      return(Edge(this->graph, i, adj_node_index, index));
    }

      /* Increments to the next edge in the current node's neighborhood
     * @return the edge iterator pointing at the next edge in the current node's neighborhood
     */
    IncidentIterator& operator++(){
      // Update the edge
      this->edge_adj_id ++;
      return(*this); 
    }


      /* Defines equality between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is equal to the argument, false otherwise
     */
    bool operator==(const IncidentIterator& it) const{
      return(this->graph == it.graph && this->node_index == it.node_index && this->edge_adj_id == it.edge_ajd_id);
    }

      /* Defines difference between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is different from the argument, false otherwise
     */
    bool operator!=(const IncidentIterator& it) const{
      return(!(this->graph == it.graph && this->node_index == it.node_index && this->edge_adj_id == it.edge_adj_id));
    }


   private:
    friend class Graph;

    // Graph in which we are
    Graph* graph;

    // Index of the current node (we iterate over its edges)
    size_type node_index;

    // Index of the current edge in the list of edges adjacency_list[node_index]
    size_type edge_adj_id;

    IncidentIterator(const Graph* g, size_type node_id, size_type edge_id){
      this->graph = const_cast<Graph*>(g);
      this->node_index = node_id;
      this->edge_adj_id = edge_id;
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
    EdgeIterator() {
    }

      /* Dereference operator
     * @return the node behind the current iterator
     */
    Edge operator*() const{
      return((*(graph)).edges[(*graph).edge_i2u[edge_index]]);
    }


      /* Increments to the next edge in the graph
     * @return the node iterator pointing at the next edge in the graph
     */
    EdgeIterator& operator++(){
      this->edge_index ++;
      return(*this);
    }

      /* Defines equality between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is equal to the argument, false otherwise
     */
    bool operator==(const EdgeIterator& it) const{
      return(this->graph == it.graph && this->edge_index == it.edge_index);
    }

      /* Defines difference between two iterators
     * @param[in, out] iterator that we compare with the current iterator
     * 
     * @return true if the current iterator is different from the argument, false otherwise
     */
    bool operator!=(const EdgeIterator& it) const{
      return(!(this->graph == it.graph && this->edge_index == it.edge_index));
    }


   private:
    friend class Graph;
    // Graph in which we are
    Graph* graph;

    // Index of the current edge in the list of edges of our graph
    size_type edge_index;

    EdgeIterator(const Graph* g, size_type edge_id){
      this->graph = const_cast<Graph*>(g);
      this->edge_index = edge_id;
    }

  };

   /* Get access to the first iterator when we iterate on the current graph's edges.
  * @return the the first iterator when we iterate on the current graph's edges
  */
  edge_iterator edge_begin() const{
    return(EdgeIterator(this, 0));
  }

   /* Get access to the last iterator when we iterate on the current graph's edges.
  * @return the the last iterator when we iterate on the current graph's edges
  */
  edge_iterator edge_end() const{
    return(EdgeIterator(this, num_edge));
  }


 private:

};

#endif // CME212_GRAPH_HPP
