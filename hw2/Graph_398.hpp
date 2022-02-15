#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>

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

  using node_value_type = V; // NodeData in HW2

  using edge_value_type = E; // EdgeData in HW2
 
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
  Graph(): nodes_(), valid_nodes_(), edges_(), valid_edges_() {
    // HW0: YOUR CODE HERE
    // using initializer list to default construct nodes_ and edges_
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
      graph_ptr = nullptr;
      idx_ = 0; 
    }

    Point& position() {
      return (*this->graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[idx_] ].node_point;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // dereference graph pointer,
      // find internal node with the node's idx
      // and return the Point object associated with it
      return (*this->graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[idx_] ].node_point; 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      //internal_node& n = (*this->graph_ptr).nodes_[this->idx_]; 
      return (*this->graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[idx_] ].node_value;
    }

    const node_value_type& value() const {
      //const internal_node& n = (*this->graph_ptr).nodes_[this->idx_]; 
      return (*this->graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[idx_] ].node_value;
    }
    
    size_type degree() const {
      return (*this->graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[idx_] ].adj_nodes_.size();
    }
    
    incident_iterator edge_begin() const {
      unsigned int begin_i = 0;
      return IncidentIterator(this->graph_ptr, this->idx_, begin_i);
    }

    incident_iterator edge_end() const{
      return IncidentIterator(this->graph_ptr, this->idx_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      return this->graph_ptr == n.graph_ptr and this->idx_ == n.idx_;
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
      //(void) n;           // Quiet compiler warning
      // assert to only compare nodes in the same graph
      //assert(this->graph_ptr == n.graph_ptr);
      // check if sum of three conditions only equate to 1
      bool trichotomy = (this->idx_ < n.idx_) + (this->idx_ == n.idx_)
	+ (this->idx_ > n.idx_);
      return (this->graph_ptr < n.graph_ptr) || 
        ( this->idx_ < n.idx_ and trichotomy == 1 ); 
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare privayte data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_ptr; // pointer to Graph class that contains Nodes
    size_type idx_; // unique identification/index for a node
    
    /*Private Constructor*/
    Node(const Graph* g, size_type idx)
        : graph_ptr(const_cast<Graph*>(g)), idx_(idx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->valid_nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& v = node_value_type() ) {
    // ASK: what are we supposed to do with node_value_type argument here?
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    
    // instantiate internal_node object
    // add to Graph's vector of internal nodes
    // return Node proxy obj that will point to the new internal node

    internal_node new_internal_node;
    new_internal_node.node_point = position;
    new_internal_node.node_value = v;
    new_internal_node.adj_nodes_ = std::vector<adjacent_nodes>();
    new_internal_node.valid_idx = valid_nodes_.size();
    this->nodes_.push_back(new_internal_node);

    valid_nodes_.push_back( nodes_.size() - 1);

    return Node(this, this->num_nodes()-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    // check  if Node's idx is smaller than the size
    // and that we are comparing the same graph 
    return n.idx_ < this->size() and this == n.graph_ptr;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    // create new Node obj with idx i;
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
      graph_ptr = nullptr;
      idx_ = 0;
      dom_idx_ = 0; 
      sub_idx_ = 0; 
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // dereference the graph,
      // get edge list, return the index
      // get index to 1st elt in a pair
      //size_type node1_idx = std::get<0>((*graph_ptr).edges_[idx_]);
      return Node(this->graph_ptr, dom_idx_);     
    }

    /* Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // derefence the graph,
      // get edge list, return the index
      // get index to 2nd elt in a pair
      //size_type node2_idx = std::get<1>((*graph_ptr).edges_[idx_]);
      // TO DO: return correct orientation of node 
      return Node(this->graph_ptr, sub_idx_);
    }

    /* Getter method to return Edge index */
    size_type get_edge_index() const {
      return this->idx_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return this->graph_ptr == e.graph_ptr and this->idx_ == e.idx_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //assert(this->graph_ptr == e.graph_ptr);
      return this->graph_ptr < e.graph_ptr || 
        this->idx_ < e.idx_;
    }

    /* Return the length of an edge */
    double length() const {
      return norm( this->node1().position() - this->node2().position() ); 
    }

    edge_value_type& value(){
      //internal_edge& e = (*this->graph_ptr).edges_[this->idx_]; 
      return (*this->graph_ptr).edges_[ (*graph_ptr).valid_edges_[idx_] ].edge_value;
    }

    const edge_value_type& value() const {
      //const internal_edge& e = (*this->graph_ptr).edges_[this->idx_]; 
      return (*this->graph_ptr).edges_[ (*graph_ptr).valid_edges_[idx_] ].edge_value;
    }
    

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_ptr; 
    size_type idx_;
    size_type dom_idx_; 
    size_type sub_idx_; 

    /* Private Constructor */
    Edge(const Graph* g, size_type idx, size_type dom, size_type sub)
        : graph_ptr(const_cast<Graph*>(g)), idx_(idx), dom_idx_(dom), sub_idx_(sub) {
    }
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->valid_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    size_type dom = this->edges_[ valid_edges_[i] ].node_idx.first;
    size_type sub = this->edges_[ valid_edges_[i] ].node_idx.second;
    return Edge(this, i, dom, sub);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning

    std::vector<adjacent_nodes> a_adj = nodes_[ valid_nodes_[a.idx_] ].adj_nodes_;
    // If edge already exists
    for(size_type i = 0; i < a.degree(); ++i){
      if( valid_nodes_[b.idx_] == a_adj[i].adj_node_idx) {
      //if(b.idx_ == a_adj[i].adj_node_idx) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e = edge_value_type()) {
    // HW0: YOUR CODE HERE
    //(void) a, (void) b;   // Quiet compiler warning
    
    // check if an Edge already exists, find the index,
    // and then return the Edge object 
    bool exists = false;
    unsigned int found_idx;

    std::vector<adjacent_nodes> a_adj = nodes_[ valid_nodes_[a.idx_] ].adj_nodes_;
    // If edge already exists
    for(unsigned int i =0; i < a_adj.size(); ++i){
      if(valid_nodes_[b.idx_] == a_adj[i].adj_node_idx ) {
        exists = true; 
        // found_idx = a_adj[i].edge_idx];
        found_idx = edges_[ a_adj[i].edge_idx ].valid_idx; 
      }
    }

    // if edge does not exist
    // create new adjacent_node objects and add to internal_node 
    if(!exists) {
      adjacent_nodes to_a, to_b; 

      to_a.edge_idx = edges_.size();
      //to_a.edge_idx = this->num_edges();
      // to_a.adj_node_idx = b.idx_; 
      to_a.adj_node_idx = b.graph_ptr->valid_nodes_[b.idx_]; 
      //a.graph_ptr->nodes_[a.idx_].adj_nodes_.push_back(to_a);
      a.graph_ptr->nodes_[ valid_nodes_[a.idx_] ].adj_nodes_.push_back(to_a);

      to_b.edge_idx = edges_.size();
      //to_b.edge_idx = this->num_edges();
      //to_b.adj_node_idx = a.idx_; 
      to_b.adj_node_idx = a.graph_ptr->valid_nodes_[a.idx_]; 
      //b.graph_ptr->nodes_[b.idx_].adj_nodes_.push_back(to_b);
      b.graph_ptr->nodes_[ valid_nodes_[b.idx_] ].adj_nodes_.push_back(to_b);

      std::pair<int, int> p(a.idx_, b.idx_);
      //edges_.push_back(p);
      internal_edge new_internal_edge;
      new_internal_edge.edge_value = e;
      new_internal_edge.node_idx = p;
      new_internal_edge.valid_idx = this->num_edges();
      // assign a new id edges_.size(); 
      edges_.push_back(new_internal_edge);

      valid_edges_.push_back( edges_.size()-1 ); 

      return Edge(this, valid_edges_.size()-1, a.idx_, b.idx_);
    } else {
      return Edge(this, found_idx, a.idx_, b.idx_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodes_.clear();
    this->edges_.clear();
    this->valid_nodes_.clear();
    this->valid_edges_.clear();
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
      graph_ptr = nullptr;
      nidx_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const{
      // NOTE: USE UID TO INDEX 
      return Node(this->graph_ptr, this->nidx_);
    }

    NodeIterator& operator++(){
      ++nidx_;
      return *this;
    }

    bool operator==(const NodeIterator& it) const {
      return this->graph_ptr == it.graph_ptr and 
        this->nidx_ == it.nidx_;
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_ptr; 
    size_type nidx_;

    /*Private Constructor*/
    NodeIterator(const Graph* g, size_type nidx)
        : graph_ptr(const_cast<Graph*>(g)), nidx_(nidx) {
    }
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    unsigned int begin_i =0;
    return NodeIterator(this, begin_i);
  }
  
  node_iterator node_end() const {
    return NodeIterator(this,this->size()); 
  }

  //
  // Incident Iterator
  //

  /** @class Graph::
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator :private totally_ordered<IncidentIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_ptr(nullptr), node1_idx_(0), adj_idx_(0) {

    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      size_type e =  (*graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[node1_idx_] ].adj_nodes_[adj_idx_].edge_idx;
      size_type a = (*graph_ptr).nodes_[ (*graph_ptr).valid_nodes_[node1_idx_] ].adj_nodes_[adj_idx_].adj_node_idx;
      return Edge(this->graph_ptr, (*graph_ptr).edges_[e].valid_idx , node1_idx_, (*graph_ptr).nodes_[a].valid_idx );
    }
    
    IncidentIterator& operator++(){
      ++adj_idx_;
      return *this;
    }

    bool operator==(const IncidentIterator& it) const {
      return this->graph_ptr == it.graph_ptr and  
        this->node1_idx_ == it.node1_idx_ and
        this->adj_idx_ == it.adj_idx_; 
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_ptr; 
    size_type node1_idx_;
    size_type adj_idx_;

    /*Private Constructor*/
    IncidentIterator(const Graph* g, size_type n_idx, size_type adj_idx)
        : graph_ptr(const_cast<Graph*>(g)), node1_idx_(n_idx), adj_idx_(adj_idx) {
    }
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_ptr = nullptr;
      eidx_ = 0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      size_type dom = (*graph_ptr).edges_[ (*graph_ptr).valid_edges_[eidx_] ].node_idx.first; 
      size_type sub = (*graph_ptr).edges_[ (*graph_ptr).valid_edges_[eidx_] ].node_idx.second; 
      return Edge(this->graph_ptr, this->eidx_, dom, sub);
    }

    EdgeIterator& operator++(){
      ++eidx_;
      return *this;
    }

    bool operator==(const EdgeIterator& it) const {
      return this->graph_ptr == it.graph_ptr and 
        this->eidx_ == it.eidx_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_ptr;
    size_type eidx_;

    /*Private Constructor*/
    EdgeIterator(const Graph* g, size_type eidx)
        : graph_ptr(const_cast<Graph*>(g)), eidx_(eidx) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const {
      unsigned int begin_i =0;
      return EdgeIterator(this, begin_i); 
   }
  
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  }

  /** Remove an edge of the graph given a (user-facing, valid) Edge index, 
   * returning 1 if successful, 0 if not.
   * @pre Edge @e e is a distinct valid edge of this graph
   * @return size_type of 1 if Edge @e e was succesfully removed, 
   * 0 if Edge @e did not exist. 
   * @post has_edge( @e e.node1(), @e e.node2()) == false
   * @post If old has_edge( @e e.node1(), @e e.node2()),
   *              new num_edges() == old num_edges() - 1.
   *       Else,  new num_edges() == old num_edges().
   * @post decreases the degree of nodes connected to @e e by 1.
   * 
   * Removes Edge @e e from the vector of valid edges, and removes 
   * corresponding adjacent node objects from e.node1() and e.node2().
   * 
   * Complexity: Max complexity is o more than O(num_edges())
   */
size_type remove_edge(const Edge & e){

  if( !has_edge(e.node1(), e.node2())) {
    return 0;
  }

  // get both user-facing (valid) index and internal unique id index
  size_type e_idx = e.get_edge_index(); 
  size_type uid = valid_edges_[e_idx];

  // remove from other nodes' (dom and sub nodes) adjacent nodes
  size_type n1_idx = e.node1().index(); 
  size_type n2_idx = e.node2().index();

  // iterate through all adjacent nodes of node1
  for(size_type i = 0; i < e.node1().degree(); ++i){
      // check if for each adjacent node, it has the edge that we want to delete 
      if (nodes_[ valid_nodes_[n1_idx] ].adj_nodes_[i].edge_idx == uid) { 
        // swap and pop
        std::iter_swap(nodes_[ valid_nodes_[n1_idx] ].adj_nodes_.begin() + i , 
          nodes_[ valid_nodes_[n1_idx] ].adj_nodes_.end() - 1 );
        nodes_[ valid_nodes_[n1_idx] ].adj_nodes_.pop_back();
    }
  }

  // iterate through all adjacent nodes of node2
  for(size_type i = 0; i < e.node2().degree(); ++i){
    // check if for each adjacent node, it has the edge that we want to delete 
    if (nodes_[ valid_nodes_[n2_idx] ].adj_nodes_[i].edge_idx == uid) { 
      // swap and pop
        std::iter_swap(nodes_[ valid_nodes_[n2_idx] ].adj_nodes_.begin() + i ,
         nodes_[ valid_nodes_[n2_idx] ].adj_nodes_.end() - 1 );
        nodes_[ valid_nodes_[n2_idx] ].adj_nodes_.pop_back();
    }
  }

  // remove from valid_edges_ with uid index 
  std::iter_swap(valid_edges_.begin() + e_idx , valid_edges_.end() - 1 );
  valid_edges_.pop_back();
  // update and re-assign the valid index
  edges_[ valid_edges_[e_idx] ].valid_idx = e_idx;

  return 1;
}

  /** Remove an edge of the graph connected by two given nodes,
   *  returning 1 if successful, 0 if not.
   * @pre Node @n1 n1 and @n2 n2 are distinct valid nodes of this graph
   * @return size_type of 1 if Edge connected by @n1 n1 and @n2 n2 was 
   * succesfully removed, 0 if such Edge did not exist. 
   * @post has_edge( @n1 n1), @n2 n2) == false
   * @post If old has_edge( @e e.node1(), @e e.node2()),
   *           new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   * @post decreases the degree of nodes connected to @e e by 1.
   * 
   * Removes all Edges connected by @n1 n1 and @n2 n2, by calling 
   * remove_edge on the Edge object iterated with the IncidentIterator, 
   * removing edges from the vector of valid edges, and removes 
   * corresponding adjacent node objects from e.node1() and e.node2().
   * 
   * Complexity: Max complexity is o more than O(num_edges())
   */
size_type remove_edge(const Node& n1, const Node & n2){
  if( !has_edge(n1, n2)) {
    return 0;
  } 

  // iterate through n1 incident nodes, check if the second node is n2, 
  // if so, remove the edge
  for (auto nit = n1.edge_begin(); nit != n1.edge_end(); ++nit ){
    if ( (*nit).node2() == n2) {
      size_type res = remove_edge(*nit); 
      return res;
    }
  }
}

  /** Remove an edge of the graph given an EdgeIterator.
   *  returning 1 if successful, 0 if not.
   * @pre EdgeIterator @e_it e_it is an iterator to valid edge of this graph
   * @return EdgeIterator of this->edge_begin()
   * @post has_edge( (*eit).node1(), (*e_it).node2() ) == false
   * @post If old has_edge( (*eit).node1(), (*e_it).node2() ),
   *           new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   * @post decreases the degree of nodes connected to @e e by 1.
   * 
   * Removes all Edges referenced by EdgeIterator @e_it e_it by calling 
   * remove_edge on the Edge object iterated with the EdgeIterator, 
   * removing edges from the vector of valid edges, and removes 
   * corresponding adjacent node objects from (*eit).node1() and (*eit).node2().
   * 
   * Complexity: Max complexity is o more than O(num_edges())
   */
edge_iterator remove_edge(edge_iterator e_it){
  remove_edge(*e_it);
  return this->edge_begin(); 
}

  /** Remove a node of the graph given a Node object.
   *  returning 1 if successful, 0 if not.
   * @pre Node @n n is a valid node of this graph
   * @return size_type of 1 if Node @n n was 
   * succesfully removed, 0 if such Node did not exist. 
   * @post has_node( @n n ) == false
   * @post If old has_node( @n n ),
   *           new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   * 
   * Removes all Edges incident to Node @n n by calling 
   * remove_edge on the Edge objects iterated with IncidentIterator.
   * Then removes Node @n n from the valid_nodes_ vector. 
   * 
   * Complexity: No more than O(num_nodes())
   */
size_type remove_node(const Node & n){
  if(!has_node(n)){
    return 0;
  }

  // while node still has incident edges, remove the incident edge
  auto iit = n.edge_begin();
  while( n.degree() > 0 ) {
    remove_edge(*iit);
  }

  // remove from valid_nodes_
  size_type n_idx = n.index(); 
  std::iter_swap(valid_nodes_.begin() + n_idx , valid_nodes_.end() - 1 );
  valid_nodes_.pop_back();
  
  // update valid_idx in nodes_ 
  nodes_[ valid_nodes_[n_idx] ].valid_idx = n_idx; 

  return 1; 
}

  /** Remove a node of the graph given a NodeIterator object.
   *  returning 1 if successful, 0 if not.
   * @pre NodeIterator @n_it n_it is a NodeIterator to a valid node of this graph
   * @return NodeIterator of this->node_begin()
   * @post has_node( *(n_it) ) == false
   * @post If old has_node( (*n_it) ),
   *           new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   * 
   * Removes Node referenced by @n_it n_it by calling 
   * remove_node on the Node object iterated with NodeIterator @n_it n_it.
   * 
   * Complexity: No more than O(num_nodes())
   */
node_iterator remove_node(node_iterator n_it){
  remove_node(*n_it);
  return this->node_begin(); 
}

 private:
 // struct to hold adjacent nodes to a particular node
  struct adjacent_nodes {
    size_type edge_idx; // edge unique id (uid)
    size_type adj_node_idx; // adjacent node unique id (uid)
  };

  // struct to hold internal node 
  struct internal_node {
    Point node_point; // holds (x,y,z) position of the node
    node_value_type node_value;
    std::vector<adjacent_nodes> adj_nodes_;
    size_type valid_idx;
  };

  // struct to hold internal edge
  struct internal_edge {
    edge_value_type edge_value;
    std::pair<int,int> node_idx;
    size_type valid_idx; 
  };

  std::vector<internal_node> nodes_;
  std::vector<size_type> valid_nodes_; 
  std::vector<internal_edge> edges_; 
  std::vector<size_type> valid_edges_; // contains edge idx 

  // HW0 YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
