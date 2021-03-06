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
template <typename V>
class Graph {
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;


  /**Type of node value */
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



 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  std::vector<std::pair<Point*, node_value_type>> nodes_;
  std::vector<std::vector<unsigned>> edges_;
  std::map<unsigned, std::vector<unsigned> > edge_map_;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //


  

  /** Construct an empty graph. */
  Graph()
    // HW0: YOUR CODE HERE
    : nodes_(std::vector<std::pair<Point*,node_value_type>>(0)),
      edges_(std::vector<std::vector<size_type>>(0)),
      edge_map_(std::map<size_type, std::vector<size_type>>()) {
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
  class Node : private totally_ordered<Node>{
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
      return *((*graph_).nodes_[node_idx_].first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(0<= node_idx_ && node_idx_ < graph_->size());
      return node_idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    /** Access this node's value
	@return this node's value of type node_value_type
     */
    node_value_type& value(){
      return (*graph_).nodes_[node_idx_].second;
    }


    /** Constant version: Access this node's value
	@return this node's value of type const node_value_type
     */
    const node_value_type& value() const {
      return (*graph_).nodes_[node_idx_].second;
    }


    /**@brief calculate the number of edges incident to the node
     * @return number of incident edges

     */
    size_type degree() const{
      return graph_->edges_[node_idx_].size();
    }

    /**@brief find the first edge related to the incident iterator
     * for this node
     *@return iterator pointing to first edge, which is the first
     *node that appears in the vector of the edges_ matrix associated
     *with this node
     */
    
    incident_iterator edge_begin() const {
      if(degree()==0){
	return edge_end();}
      
      incident_iterator it = IncidentIterator(graph_, node_idx_,0);
      return it;
    }

    /**@brief find the last edge related to the incident iterator
     * for this node
     *@return iterator pointing to end the iteration
     */
    
    incident_iterator edge_end() const {
      incident_iterator it = IncidentIterator(graph_,node_idx_,degree());
      return it;
    }
    
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(n.graph_ == graph_ && n.node_idx_ == node_idx_){
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
      // HW0: YOUR CODE HERE
      if(node_idx_< n.node_idx_) {
	return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph 
    Graph* graph_;
    // This node's index number
    size_type node_idx_;
    Node(const Graph* graph, size_type node_idx)
      : graph_(const_cast<Graph*>(graph)), node_idx_(node_idx) {
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
   *No value is specified for this node
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    Point* newpoint = new Point;
    *newpoint = position;
    node_value_type new_val = node_value_type();
    std::pair<Point*, node_value_type> next_node(newpoint, new_val);
    nodes_.push_back(next_node);

    //need to add new node to adjacency matrix
    std::vector<size_type> new_node;
    edges_.push_back(new_node);
   
    return Node(this, size()-1);   
  }



  /** Add a node to the graph, returning the added node.
   * A value is specified for this node
   * @param[in] position The new node's position
   * @param[in] new_val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() = new_value
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& new_val) {
    // HW0: YOUR CODE HERE
    Point* newpoint = new Point;
    *newpoint = position;
    std::pair<Point*, node_value_type> next_node(newpoint, new_val);
    nodes_.push_back(next_node);

    //need to add new node to adjacency matrix
    std::vector<size_type> new_node;
    edges_.push_back(new_node);
   
    return Node(this, size()-1);   
  }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(this  == n.graph_){
      return true;
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
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
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
  class Edge : private totally_ordered <Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_idx_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_idx_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if((e.node1_idx_ == node1_idx_ && e.node2_idx_ == node2_idx_)||
	 (e.node1_idx_ == node2_idx_ &&
	  e.node2_idx_ == node1_idx_)) return true;     
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      size_type min_e = std::min(e.node1_idx_, e.node2_idx_);
      size_type max_e = std::max(e.node1_idx_, e.node2_idx_);
      size_type min_this = std::min(node1_idx_, node2_idx_);
      size_type max_this = std::max(node1_idx_, node2_idx_);
      if((min_this<min_e) || ((min_this == min_e)
			     && (max_this < max_e)) ) return true;
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
    size_type node1_idx_;
    size_type node2_idx_;
    
    Edge(const Graph* graph, size_type node1_idx, size_type node2_idx)
      : graph_(const_cast<Graph*>(graph)), node1_idx_(node1_idx),
	node2_idx_(node2_idx){
    }
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    
    return edge_map_.size();//num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */


  
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    
    return Edge(this, edge_map_.at(i)[0], edge_map_.at(i)[1]);
  }
  


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(this == a.graph_ && this == b.graph_);
    size_type a_idx = a.node_idx_;
    size_type b_idx = b.node_idx_;
    assert(a_idx < size() && b_idx < size());

    for (auto & element : edges_[a_idx]){
      if(element == b_idx){
	return true;
      }
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
    assert(this == a.graph_ && this == b.graph_);
    size_type a_idx = a.node_idx_;
    size_type b_idx = b.node_idx_;
    assert(a_idx < size() && b_idx < size() && a_idx != b_idx);

    
    for (auto & element : edges_[a_idx]){
      if(element == b_idx){
	return Edge(this, a_idx, b_idx);
      }
    }
    
    
 
    edges_[a_idx].push_back(b_idx);
    edges_[b_idx].push_back(a_idx);

    std::vector<size_type> edge_vec(2);
    edge_vec[0] = a_idx;
    edge_vec[1] = b_idx;
    size_type eid = this->num_edges();
    edge_map_[eid] = edge_vec;

    
    //assert(this->has_edge(a,b) == true && has_edge(b,a) == true);
    //assert(this->num_edges() == eid+1);
    return Edge(this, a_idx, b_idx);  
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
    edge_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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


    /**@brief Define dereferencing operator
     * @return the node being referenced*/
    Node operator*() const{
      return Node(graph_, cur_node_);
    }


    /**@brief Define increment  operator
     * @return the node iterator
     * @post cur_node_ 1 greater than before operator
     * was applied */
    node_iterator& operator++(){
      ++ cur_node_;
      return *this;
    }


    /**@brief Define equality  operator
     *@param[in] it iterator
     * @return whether or not the iterators are equal*/
    bool operator==(const NodeIterator& it) const{
      return graph_ == it.graph_ && cur_node_ == it.cur_node_;
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type cur_node_;
    NodeIterator(const Graph* graph, size_type cur_node) :
      graph_(const_cast<Graph*>(graph)), cur_node_(cur_node){}
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /**@brief Define begin function for node iterator 
     * @return iterator that points to beginning of nodes container*/
  node_iterator node_begin() const {
    node_iterator it = NodeIterator(this,0);
    return it;
  }

  /**@brief Define end function for iterator 
     * @return iterator that points to end of nodes container*/
  node_iterator node_end() const {
    node_iterator it = NodeIterator(this, size());
    return it;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator> {
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


    /**@brief Define dereferencing operator
     * @return the edge being referenced*/
    Edge operator*()const{
      return Edge(graph_, start_node_, graph_->edges_[start_node_][inc_idx_]); 
    }


    /**@brief Define increment  operator
     * @return the incident iterator
     * @post inc_node_ is increased by 1 and start_node_ stays same */
    IncidentIterator& operator++(){
      ++inc_idx_;
      return *this;
    }
 

    /**@brief Define equality  operator
     *@param[in] it iterator
     * @return whether or not the iterators are equal*/
    bool operator ==(const IncidentIterator& it) const {
      return (graph_ == it.graph_ && start_node_ == it.start_node_
	      && inc_idx_ == it.inc_idx_);
    }

    
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /** Keep track of the node with which we are iterating over
	incident edges, and the corresponding node for the current
	edge
     */
    Graph* graph_;
    size_type start_node_;
    size_type inc_idx_;
    IncidentIterator(const Graph* graph,
		     size_type start_node, size_type inc_idx) :
      graph_(const_cast<Graph*>(graph)), start_node_(start_node),inc_idx_(inc_idx){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /**@brief define edge iterator dereferencing operator
     @return current edge
     */
    Edge operator*() const{
      size_type a_idx = graph_ ->edge_map_.at(edge_id_)[0];
      size_type b_idx = graph_ ->edge_map_.at(edge_id_)[1];
      return Edge(graph_, a_idx, b_idx);
    }
    
    /**@brief define edge iterator increment operator
     *@return the iterator with the id incremented
     */
    
    EdgeIterator& operator++(){
      ++edge_id_;
      return *this;
    }

    /**@brief define equality operator for edge iterator 
     * @param[in] an edge iterator to be compared to
     */
    bool operator==(const EdgeIterator& edgeit) const{
      return (graph_ == edgeit.graph_ && edge_id_ == edgeit.edge_id_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_id_;

    /**Construct a valid EdgeIterator */
    EdgeIterator(const Graph* graph, size_type edge_id)
      : graph_(const_cast<Graph*>(graph)), edge_id_(edge_id){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:


  /**@brief find the first edge to begin the iteration
   *@return the first edge id (0) which serves as a key to 
   *the edge_map_
   */
  edge_iterator edge_begin() const{
    edge_iterator it = EdgeIterator(this, 0);
    return it;
  }
  

  /**@brief supply the ending condition to the iterator
   *@return the size of the edge_map_ 
   */
  edge_iterator edge_end() const{
    edge_iterator it = EdgeIterator(this, num_edges());
    return it;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP

