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

#include <unordered_map>


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> 
class Graph {

 // Predeclare the internal structs
 struct internal_node;
 struct internal_edge; 

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  std::vector<internal_node> nodes_; // a vector of internal nodes 
  std::vector<internal_edge> edges_; // a vector of internal edges 

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

  using node_value_type = V;
  using edge_value_type = E; 

  using index_umap = std::unordered_map<size_type, size_type>; 


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  // HW0: YOUR CODE HERE
  Graph() 
    : nodes_() {
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
     // HW0: YOUR CODE HERE
    Node() {
    }

    /** Return this node's position (to be modified later) */
    Point& position(){
      return fetch().point;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().uid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Returns the value of the node 
    * @return The value of the current node 
    */
    node_value_type& value(){
      return fetch().value; 
    };

    /** Returns the value of a node (as const) 
    * @return The value of the current node 
    */
    const node_value_type& value() const{
      return fetch().value; 
    };

    /** Returns the degree of this node 
    * @return The degree of the current node 
    */
    size_type degree() const{
      return fetch().incident_edges.size();
    };

    /** Gets the beginning of the incident edge iterator
    * @return An IncidentIterator pointing to the first incident edge of this node
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, this, fetch().incident_edges.begin());
    };
    
    /** Gets the end of the incident edge iterator 
    * @return An IncidentIterator pointing to one past the last incident edge of this node
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, this, fetch().incident_edges.end());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      
      bool same_uid = (uid_ == n.uid_); 
      bool same_graph = (graph_ == n.graph_);

      return same_uid && same_graph;
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
      (void) n;           // Quiet compiler warning
      return uid_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container 
    Graph* graph_;
    // This node's unique identification number
    size_type uid_;

    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    /** Helper method to return the appropriate node.
     * Returns the node with the correct uid
     */
    internal_node& fetch() const {
      return graph_->nodes_[uid_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning

    // the uid for our next node 
    size_type next_uid = num_nodes();

    // create a new internal node 
    internal_node next_internal_node; 
    next_internal_node.point = position; 
    next_internal_node.uid = next_uid; 
    next_internal_node.value = value;
    nodes_.push_back(next_internal_node);

    // create the new node 
    return Node(this, next_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    
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
    // HW0: YOUR CODE HERE
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return fetch_edge().node1; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch_edge().node2;  
    }

    /** Return the length of this edge */
    double length() const{
      return graph_->euclidean_dist_(node1(), node2()); 
    }

    /** Return the value of this edge */ 
    edge_value_type& value(){
      return fetch_edge().value;
    }

    /** Return the value of this edge */ 
    const edge_value_type& value() const{
      return fetch_edge().value; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE

      // test one direction 
      bool direction1 = (node1() == e.node1()) && (node2() == e.node2()); 
      // test other direction 
      bool direction2 = (node1() == e.node2()) && (node2() == e.node1()); 

      // return if either is true 
      return (direction1 || direction2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return node1().index() < e.node1().index(); 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_; 
    size_type uid_; 

    Edge(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }    

    /** Helper method to return the appropriate edge.
     * Returns the edge with the correct uid
     */
    internal_edge& fetch_edge() const {
      return graph_->edges_[uid_];
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
    (void) i;             // Quiet compiler warning
    return Edge(this, i); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (incident_iterator it = a.edge_begin(); it != a.edge_end(); ++it){
      if((*it).node2() == b){
        return true; 
      }
    }
    return false; 
  }

  // like has_edge(), except it returns the index of the Edge it finds
  // it returns -1 if it can't find the edge 
  int find_edge(const Node& a, const Node& b) const {
    for (incident_iterator it = a.edge_begin(); it != a.edge_end(); ++it){
      if((*it).node2() == b){
        return (int) (*it).uid_; 
      }
    }
    return -1; 
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    (void) a, (void) b;   // Quiet compiler warning

    // see if this edge is already in the graph 
    int edge_index = find_edge(a, b);
    if (edge_index != -1){
      return Edge(this, edge_index);
    }

    // otherwise we create a new edge with these nodes 
    size_type next_uid = num_edges();

    // create a new internal edge 
    internal_edge next_internal_edge; 
    next_internal_edge.node1 = a;   
    next_internal_edge.node2 = b; 
    next_internal_edge.uid = next_uid; 
    next_internal_edge.value = value; 
    edges_.push_back(next_internal_edge);

    // update the internal nodes that are a part of this edge
    // the key is the uid of the node, the value is it's index in edges_
    nodes_[a.index()].incident_edges[next_uid] = next_uid;     
    nodes_[b.index()].incident_edges[next_uid] = next_uid;   

    // create the new edge 
    return Edge(this, next_uid);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{

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

    /** Increments the NodeIterator 
    * @return The next node for this iterator 
    * @post index_ is increased by 1 
    */
    NodeIterator& operator++(){
      index_++; 
      return *this; 
    }    

    /** Compares two NodeIterator objects
    * @param[in] other_node NodeIterator 
    * @return true if this == other_node_iter 
    */
    bool operator==(const NodeIterator& other_node_iter) const{
      bool same_graph = (graph_ == other_node_iter.graph_);
      bool same_index = (index_ == other_node_iter.index_);
      return same_graph && same_index; 
    }

    /** Dereferences the node that NodeIterator points to  
    * @return The node this iterator points to 
    */
    Node operator*() const{
      return Node(graph_, index_);
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // pointer to the Graph 
    size_type index_; // index of the current node 

    //Private constructor that can be accessed by Graph 
        NodeIterator(const Graph* graph, size_type index) 
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a NodeIterator pointing to the first Node in Graph
  * @return NodeIterator pointing to Node 0 
  */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Returns a NodeIterator pointing to one place after the last Node in Graph
  * @return NodeIterator point to one past the last Node 
  */
  node_iterator node_end() const{
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    /** Increments the IncidentIterator
    * @return This IncidentIterator pointing to the next incident edge
    * @post index_ is increased by 1 
    */
    IncidentIterator& operator++(){

      it_++; 
      return *this; 

    };

    /** Checks if two IncidentIterators are the same 
    * @param[in] other_incident_iterator Another IncidentIterator
    * @return true if this == other_incident_iterator
    */ 
    bool operator==(const IncidentIterator& other_incident_iterator) const{

      bool same_graph = (graph_ == other_incident_iterator.graph_);
      bool same_node = (node_ == other_incident_iterator.node_); 
      bool same_it = (it_ == other_incident_iterator.it_); 

      return (same_graph && same_node && same_it); 

    };

    /** Returns the current Edge in this IncidentIterator 
    * @return The Edge this IncidentIterator points to 
    * @post For the returned edge, e, it will always hold that e.node1()==node_ 
    */
    Edge operator*() const{

      // get the edge in edges_ 
      int index = it_->second; 
      internal_edge this_edge = graph_->edges_[index];

      // if the first node in this edge is the node that is calling the IncidentIterator
      // then we return this edge, as it is already "oriented"
      if (this_edge.node1.index() == node_->index()){
        return Edge(graph_, index);
      }; 

      // otherwise, we need to manually flip the edge so it's properly oriented 
      internal_edge flipped_edge = this_edge; 
      flipped_edge.node1 = this_edge.node2; 
      flipped_edge.node2 = this_edge.node1; 
      graph_->edges_[index] = flipped_edge;

      return Edge(graph_, index);
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    Graph* graph_; // pointer to the Graph 
    Node* node_; // pointer to the Node 
    index_umap::iterator it_; // umap iterator

    //Private constructor that can be accessed by the Node class.
    IncidentIterator(const Graph* graph, const Node* node, index_umap::iterator it) 
        : graph_(const_cast<Graph*>(graph)), node_(const_cast<Node*>(node)), it_(it) {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Increments the EdgeIterator
    */
    EdgeIterator& operator++(){
      index_ ++; 
      return *this; 
    };

    /** Checks if two EdgeIterators are equal 
    */ 
    bool operator==(const EdgeIterator& other_edge) const{
      bool same_graph = (graph_ == other_edge.graph_);
      bool same_index = (index_ == other_edge.index_);
      return same_graph && same_index; 
    };

    /** Dereferences this EdgeIterator 
    */ 
    Edge operator*() const{
      return Edge(graph_, index_);
    };

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    Graph* graph_; // pointer to the Graph 
    size_type index_; // index of the current edge 

    //Private constructor that can be accessed by Graph 
      EdgeIterator(const Graph* graph, size_type index) 
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns the start of the EdgeIterator
  * @return An EdgeIterator pointing to the first incident Edge 
  */ 
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0); 
  };

  /** Returns the end of an EdgeIterator
  * @return An EdgeIterator pointing to one position past the last incident Edge 
  */ 
  edge_iterator edge_end() const{
    return EdgeIterator(this, edges_.size()); 
  };

  /** Remove a node from the graph.
   * @pre @a a is a valid node of this graph
   * @return 1 if the node was removed, 0 otherwise 
   * @post has_node(@a a) == false
   * Complexity: No more than O(num nodes())
   */
  size_type remove_node(const Node& a){

  };

  /** Remove a node from the graph.
   * @pre @a n_it is a node iterator pointing to the node to be removed 
   * @return a node iterator pointing to the beginning of the graph 
   * @post has_node(*n_it) == false
   * Complexity: No more than O(num nodes())
   */
  node_iterator remove_node(node_iterator n_it){

  };

  /** Remove a node from the graph.
   * @pre @a a and @a b are valid nodes of this graph
   * @return 1 if the edge was removed, 0 otherwise 
   * @post has_edge(a, b) == false
   * Complexity: No more than O(num nodes() + num edges()) 
   */
  size_type remove_edge(const Node& a, const Node& b){
    
    // if the edge doesn't exist, return 0 
    if (!has_edge(a, b)){
      return 0; 
    }

    size_type edge_index; 
    size_type edge_uid; 
    size_type counter = 0; 

    // find the index in edges_ with our edge 
    for (internal_edge e : edges_){
      bool one_direction = (e.node1==a) && (e.node2==b); 
      bool other_direction = (e.node1==b) && (e.node2==a); 
      if(one_direction || other_direction){
        edge_index = counter;
        edge_uid = e.uid; 
        std::swap(edges_.back(), edges_[counter]);
        edges_.pop_back();
        break;
      }
      counter++; 
    }
    
    // remove from incident edge for node a
    index_umap& incident_a = nodes_[a.index()].incident_edges; 
    incident_a.erase(edge_uid); 

    for(auto it=incident_a.begin(); it!=incident_a.end(); ++it){
      if(it->second > edge_index)
        it->second--;  
    }

    // remove from incident edge for node b
    index_umap& incident_b = nodes_[b.index()].incident_edges; 
    incident_b.erase(edge_uid); 

    for(auto it=incident_b.begin(); it!=incident_b.end(); ++it){
      if(it->second > edge_index)
        it->second--;  
    }

    return 1; 

  };

  /** Remove a node from the graph.
   * @pre @a e is a valid edge of this graph
   * @return 1 if the edge was removed, 0 otherwise 
   * @post has_edge(e.node1(), e.node2()) == false
   * Complexity: No more than O(num nodes() + num edges()) 
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  };

  /** Remove a node from the graph.
   * @pre @a e_it is an edge iterator pointing to the edge to be removed
   * @return an edge iterator pointing to the beginning of the graph  
   * @post has_edge((*e_it).node1(), (*e_it).node2()) == false
   * Complexity: No more than O(num nodes() + num edges()) 
   */
  edge_iterator remove_edge(edge_iterator e_it){
    
    Edge e = *e_it; 
    size_type removed = remove_edge(e.node1(), e.node2());
    return edge_begin();
  };

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal type for nodes 
  struct internal_node {
    Point point;   // The point represented by that node 
    size_type uid; // The unique identification for that node 
    node_value_type value; // The value of this node 
    index_umap incident_edges; // The edges incident to this node 
  };

  // Internal type for edges 
  struct internal_edge {
    Node node1; // The first node of this edge
    Node node2; // The second node of this edge 
    size_type uid; // The unique identification for this edge 
    edge_value_type value; // The value of this edge 
  };

  // Calculate the euclidean distance between two nodes
  double euclidean_dist_(const Node& a, const Node& b){
    return norm(a.position() - b.position());
  }

};

#endif // CME212_GRAPH_HPP
