#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include <set>
#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <unordered_set>

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
 //Declare internal structs
 struct internal_node;
 struct internal_edge;
 struct node_adjacent;

 private:
  // HW0: CODE ADDED HERE
  // Use this space for declarations of important internal types you need
  // Interal type for nodes:
     struct internal_node{
	 Point intpoint_; //Point of node
         unsigned intnidx_; //Inward-facing Index of node
         unsigned onidx_; //Outward-facing Index of node
         V intval_; //"Value" of node
         std::vector<node_adjacent> intadjnodes_; //Vector of adjacent nodes
         internal_node(Point intpoint, unsigned intnidx, unsigned onidx, V intval) 
               : intpoint_(intpoint), intnidx_(intnidx), onidx_(onidx), intval_(intval) {}
      };
     unsigned size_; //Number of nodes
     // Internal type of edges:
     struct internal_edge{
         unsigned inteidx_; //Inward-facing Index of edge
         unsigned oeidx_; //Outward-facing Index of edge
         unsigned intnode1idx_; //Index of node1
         unsigned intnode2idx_; //Index of node2
         E inteval_;
         internal_edge(unsigned inteidx, unsigned oeidx, unsigned intnode1idx, unsigned intnode2idx, E inteval)
         :inteidx_(inteidx), oeidx_(oeidx), intnode1idx_(intnode1idx),intnode2idx_(intnode2idx),inteval_(inteval){}
     };
     unsigned num_edges_; //Number of edges
     std::vector<internal_edge> edges_vec; //Vector of internal edge objects
     std::vector<unsigned> ni2u_; //Vector of node outward-to-unique indices
     std::vector<internal_node> nodes_vec; //Vector of internal node objects
     std::vector<unsigned> ei2u_; //Vector of edge outward-to-unique indices
     //Internal type of adjacent node
     struct node_adjacent{
         unsigned adjeidx_; //Index of edge
         unsigned adjnodeidx_; //Index of node attached
         node_adjacent(unsigned adjeidx, unsigned adjnodeidx)
             : adjeidx_(adjeidx), adjnodeidx_(adjnodeidx) {}
     };
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
//
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  using node_value_type = V;
  using edge_value_type = E;
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
  Graph() 
      : size_(0), num_edges_(0) {
    // HW0: CODE ADDED HERE
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
    Node() 
        : ngraph_(nullptr), nidx_(0){
      // HW0: CODE ADDED HERE
    }
    /** Return this node's position. */
    const Point& position() const {
      // HW0: CODE ADDED HERE
      // Reach into vector of internal nodes and return point attribute.
        return ngraph_ -> nodes_vec[nidx_].intpoint_;
    }
    Point& position() {
        return ngraph_->nodes_vec[nidx_].intpoint_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    unsigned index() const {
      // HW0: CODE ADDED HERE
      // Return index attribute
      //  return nidx_;
      return ngraph_->nodes_vec[nidx_].onidx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
        return ngraph_-> nodes_vec[nidx_].intval_; //Extract the intval_
    }
    const node_value_type& value() const{
        return ngraph_-> nodes_vec[nidx_].intval_; //Extract the intval_ 
    }
    unsigned degree() const {
        // Want number of valid adjacent objects
        internal_node intnode = ngraph_-> nodes_vec[nidx_];
        return intnode.intadjnodes_.size(); //Extract number of adjacent nodes 
    }
    incident_iterator edge_begin() const {
        // Returns incident iterator pointing at front of adjacent node vector
        return IncidentIterator(ngraph_, nidx_, 0);
    }
    incident_iterator edge_end() const {
        // Returns incident iterator pointing to back of adjacent node vector
        unsigned deg = this->degree();
        return IncidentIterator(ngraph_, nidx_, deg);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: CODE ADDED HERE
      // Given 2 Nodes, return True if they have same Graph and index.
      if(this-> ngraph_ == n.ngraph_ && this-> nidx_ == n.nidx_){
          return true;
      }
      (void) n;          // Quiet compiler warning
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
      // HW0: CODE ADDED HERE
      // Overload operator to compare node indices
        //if(ngraph_ != n.ngraph_){return true;} //if graphs don't match, arbitrarily return true 
        return (ngraph_ < n.ngraph_ || this-> nidx_ < n.nidx_); 
        //(void) n;           // Quiet compiler warning
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: CODE ADDED HERE
    Graph* ngraph_; //Pointer to graph
    unsigned nidx_; //Index of node
    // Private constructor to return a valid Node
    Node(const Graph* ngraph, size_type nidx)
        : ngraph_(const_cast<Graph*>(ngraph)), nidx_(nidx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  unsigned size() const {
    // HW0: CODE ADDED HERE
    // Return size_ (number of nodes) attribute
    // return size_;
    return ni2u_.size();
    // return 0;
  }

  /** Synonym for size(). */
  unsigned num_nodes() const {
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
    // HW0: CODE ADDED HERE
    // Add a new instance of internal_node to the back of the vector
    unsigned uid = nodes_vec.size(); //Unique ID IS THIS RIGHT??? PLEASE CHECK
    unsigned oid = ni2u_.size(); //Place in active node list
    nodes_vec.emplace_back(internal_node(position, uid, oid, value)); //modified to add value 
    // Record new uid to back of ni2u_
    ni2u_.push_back(uid);
    size_ = size_ + 1; // Increment number of nodes
    //(void) position;      // Quiet compiler warning
    return Node(this, nodes_vec.size()-1);        // Return node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: CODE ADDED HERE
    // Check if node is present in the vector of nodes
    //unsigned n_idx = n.index();
    unsigned uid = ni2u_[n.index()]; 
    // See if node uid is present in active nodes
    for(unsigned i=0; i<ni2u_.size();i++){
        if(ni2u_[i] == uid){
            return true;
        }
    }
    //(void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: CODE ADDED HERE
    //(void) i;             // Quiet compiler warning
    return Node(this, ni2u_[i]);        // Return node at index i
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
    Edge()
        : egraph_(nullptr), eidx_(0), n1idx_(0), n2idx_(0) {
      // HW0: CODE ADDED HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: CODE ADDED HERE
      return Node(egraph_, n1idx_); //Return the first Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: CODE ADDED HERE
      return Node(egraph_, n2idx_); // Return the second Node
    }

    /** Return the edge length of this Edge */
    double length() const {
      Node n1 = this->node1();
      Node n2 = this->node2();
      //Get the positions of the two connected nodes
      Point pos1 = n1.position();
      Point pos2 = n2.position();
      //Compute the distance between the two positions.
      double dist = std::sqrt(std::pow(pos1.x-pos2.x,2) + 
                    std::pow(pos1.y-pos2.y,2) +
                    std::pow(pos1.z-pos2.z,2)); 
      return dist;    
    }

    /** Return value of this Edge */
    edge_value_type& value() {
        //Get the value associated with Edge
        return egraph_->edges_vec[eidx_].inteval_;
    }

    const edge_value_type& value() const {
        //Get the value associated with Edge
        return egraph_->edges_vec[eidx_].inteval_;
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: CODE ADDED HERE
      // If the graph and index of the edge match, they are ==
      if (egraph_ == e.egraph_ && eidx_ == e.eidx_){
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
      (void) e;           // Quiet compiler warning
      //HW0: CODE ADDED HERE
      //Overload operator to compare edge indices
      //if (egraph_ != e.egraph_){return true;} //If not from same graph...arbitrarily return true?
      return (egraph_ < e.egraph_ || this-> eidx_ < e.eidx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: CODE ADDED HERE
    Graph* egraph_; // Pointer to graph
    unsigned eidx_; // Index of edge
    unsigned n1idx_; // Index of first node
    unsigned n2idx_; //Index of second node
    //Private constructor for returning valid Edge objects:
    Edge(const Graph* graph, unsigned eidx, unsigned n1idx, unsigned n2idx)
    :egraph_(const_cast<Graph*>(graph)),eidx_(eidx),n1idx_(n1idx),n2idx_(n2idx){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: CODE ADDED HERE
    //  return num_edges_; //Return number of edges
    return ei2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: CODE ADDED HERE
    //(void) i;             // Quiet compiler warning
    unsigned node1idx = edges_vec[ei2u_[i]].intnode1idx_;
    unsigned node2idx = edges_vec[ei2u_[i]].intnode2idx_;
    return Edge(this, ei2u_[i], node1idx, node2idx); // Return Edge at i
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: CODE ADDED HERE
    // For each node in the corresponding adjacent nodes vector,
      // check to see if the index matches node b
    unsigned a_uid = ni2u_[a.index()];
    unsigned b_uid = ni2u_[b.index()];
    // We loop again to get the index
    for(unsigned i=0; i< nodes_vec[a_uid].intadjnodes_.size();i++){
        if(nodes_vec[a_uid].intadjnodes_[i].adjnodeidx_ == b_uid){
            return true;   
    	} 
    }  
    //(void) a; (void) b;   // Quiet compiler warning
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Check if edge already exists
    unsigned a_uid = ni2u_[a.index()];
    unsigned b_uid = ni2u_[b.index()];
    // We loop again (instead of calling has_edge) to get the index
    for(unsigned i=0; i< nodes_vec[a_uid].intadjnodes_.size();i++){
        if(nodes_vec[a_uid].intadjnodes_[i].adjnodeidx_ == b_uid){
            return Edge(this, nodes_vec[a_uid].intadjnodes_[i].adjeidx_, a_uid, b_uid);  
        }
    }     
    // If the edge doesn't exist, add a new edge
    internal_node inode_a = nodes_vec[a_uid];
    internal_node inode_b = nodes_vec[b_uid];
    unsigned e_uid = edges_vec.size();
    unsigned e_oid = ei2u_.size();

    // Add two new instances of adjacent nodes (one for a, one for b)
    // Add these instances to the corresponding vector
    nodes_vec[a_uid].intadjnodes_.push_back(node_adjacent(e_uid,b_uid));
    nodes_vec[b_uid].intadjnodes_.push_back(node_adjacent(e_uid,a_uid));
    
    // Add one new instance of internal edge and add to vector
    edges_vec.push_back(internal_edge(e_uid, e_oid, a_uid, b_uid, value));
    ei2u_.push_back(e_uid); 
   // Increment number of edges
    num_edges_ = num_edges_ + 1;
    //(void) a, (void) b;   // Quiet compiler warning
    return Edge(this, edges_vec.size()-1, a_uid, b_uid); // Return New Edge
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: CODE ADDED HERE
     ei2u_.clear();
     ni2u_.clear();
     nodes_vec.clear();
     edges_vec.clear();
     size_ = 0;
     num_edges_ = 0;
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
    //NOTE: CHANGED TO FORWARD FROM INPUT
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
        // Return the Node we are at
        internal_node intnode = it_graph_->nodes_vec.at(it_graph_->ni2u_[it_idx_]);
        unsigned uid = intnode.intnidx_; 
        return Node(it_graph_,uid); 
    }
    NodeIterator& operator++(){
        // Increment iterator index and return instance      
        it_idx_ = it_idx_ + 1;
        return *this; 
    }
    bool operator==(const NodeIterator& iterator2) const{
        // Test if two node iterators are the same
        if ((it_graph_==iterator2.it_graph_) 
          && 
         (it_idx_ == iterator2.it_idx_)){
            return true;
        }
        return false; 
   }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Pointer to graph;
    Graph* it_graph_;
    // Iterator pointer:
    unsigned it_idx_;
    // Pointer to vector of nodes
    //const std::vector<internal_node>* nvec_ptr_;
    // Private constructor:
    NodeIterator(const Graph* it_graph, unsigned it_idx)
        : it_graph_(const_cast<Graph*>(it_graph)), it_idx_(it_idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    // Return iterator pointing to front of nodes vector
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    // Return iterator pointing to back of nodes vector
    return NodeIterator(this,  ni2u_.size());
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
    //NOTE: CHANGED TO FORWARD FROM INPUT
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
        // Return the Edge the iterator is at
        // Reach into adjacency vector for information on IDs
        return Edge(graph_, 
         graph_-> nodes_vec[node1idx_].intadjnodes_[iit_idx_].adjeidx_, 
         node1idx_, 
         graph_-> nodes_vec[node1idx_].intadjnodes_[iit_idx_].adjnodeidx_);
    }
    IncidentIterator& operator++(){
        // Increment index of iterator
        iit_idx_ = iit_idx_ + 1;
        return *this;
    }
    bool operator==(const IncidentIterator& iterator2) const {
        // Test if two iterators are the same
        if((graph_==iterator2.graph_) &&
         (node1idx_==iterator2.node1idx_) &&
         (iit_idx_ == iterator2.iit_idx_)){
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Pointer to graph
    Graph* graph_;
    // Index of node of interest
    unsigned node1idx_;
    // Index of the iterator
    unsigned iit_idx_;
    // Private constructor:
    IncidentIterator(const Graph* graph, unsigned node1idx, unsigned iit_idx)
        : graph_(const_cast<Graph*>(graph)), node1idx_(node1idx), iit_idx_(iit_idx) {}
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
    //NOTE: CHANGED TO FORWARD FROM INPUT
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
        // Return the edge our iterator is at
        internal_edge intedge = eit_graph_->edges_vec.at(eit_graph_->ei2u_[eit_idx_]);
        unsigned uid = intedge.inteidx_; //Get edge index
        return Edge(eit_graph_, uid, intedge.intnode1idx_, intedge.intnode2idx_);
    }
    EdgeIterator& operator++() {
        // Increment index and return the instance
        eit_idx_ = eit_idx_ + 1;
        return *this;
    }
    bool operator==(const EdgeIterator& iterator2) const {
        // Check if two iterators are the same
        if ((eit_graph_== iterator2.eit_graph_) &&
         (eit_idx_==iterator2.eit_idx_)) {
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer to graph
    Graph* eit_graph_;
    // Index of the iterator
    unsigned eit_idx_;
    // Private constructor:
    EdgeIterator(const Graph* eit_graph, unsigned eit_idx)
    :eit_graph_(const_cast<Graph*>(eit_graph)),eit_idx_(eit_idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const {
       // Return iterator instance pointing to front of edge vector
       return EdgeIterator(this, 0);
   }
   edge_iterator edge_end() const {
       // Return iterator instance pointing to back of edge vector
       return EdgeIterator(this, ei2u_.size());
   }  

  /** @brief Remove edge @a e from the graph. 
    * @pre Edge @a e is a valid edge of the graph.
    * @return 1 to signal successful removal.
    * @post Edge @a e is removed from the graph. 
    * @post Nodes touching @a e are erased as neighbors.
    * @post Number of active edges is reduced by 1. 
    * Complexity : no more than O(num_edges() + num_nodes()) */
  size_type remove_edge(const Edge& e){
      unsigned uid = e.eidx_;
      unsigned oid = edges_vec[uid].oeidx_;
      //Erase neighbors
      unsigned n1_uid = ni2u_[e.node1().index()];
      unsigned n2_uid = ni2u_[e.node2().index()];
      for (unsigned i = 0; i < nodes_vec[n1_uid].intadjnodes_.size() ; i++){
          if (nodes_vec[n1_uid].intadjnodes_[i].adjnodeidx_ == n2_uid){
              //Erase entry
              nodes_vec[n1_uid].intadjnodes_.erase(nodes_vec[n1_uid].intadjnodes_.begin()+i);
          }
      }
      //Do the same for the other node
      for (unsigned i = 0; i < nodes_vec[n2_uid].intadjnodes_.size() ; i++){
          if (nodes_vec[n2_uid].intadjnodes_[i].adjnodeidx_ == n1_uid){
              //Erase entry
              nodes_vec[n2_uid].intadjnodes_.erase(nodes_vec[n2_uid].intadjnodes_.begin()+i);
          }
      }
      //Swap and pop
      std::iter_swap(ei2u_.begin()+oid, ei2u_.end()-1);
      ei2u_.pop_back();
      // Set last edge's new oid
      //edges_vec[ei2u_[ei2u_.size()-1]].oeidx_ = oid;
      edges_vec[ei2u_[oid]].oeidx_ = oid;
      num_edges_ = num_edges_-1;
      return 1;
  }

  /** Remove edge connecting @a n1 and @a n2 from the graph.
   * @brief Finds the edge connecting the nodes and calls remove_edge.
   * @pre Nodes @a n1 and @a n2 are valid nodes of the graph.
   * @pre Nodes @a n1 and @a n2 share an edge.
   * @return 1 if edge is successfully removed
   * @post Edge is removed from the graph.
   * @post Nodes @a n1 and @a n2 are erased as neighbors.
   * @post Number of active edges is reduced by 1. 
   * Complexity : no more than O(num_edges() + num_nodes()) */
  size_type remove_edge(const Node& n1, const Node& n2) {
     unsigned n1_uid = ni2u_[n1.index()];
     unsigned n2_uid = ni2u_[n2.index()];
     std::vector<node_adjacent> adj = nodes_vec[n1_uid].intadjnodes_;
     for(unsigned i = 0; i < adj.size(); i++){
         if (adj[i].adjnodeidx_==n2_uid){
             unsigned e_uid = adj[i].adjeidx_;
             return remove_edge(Edge(this, e_uid, n1_uid, n2_uid));
         }
     } 
     return 0;
  }

  /** Remove edge pointed at by @a e_it.
    * @brief Deferences and calls remove_edge on resulting edge.
    * @pre Edge is a valid edge of the graph.
    * @return an iterator pointing to the first edge. 
    * @post Edge is removed from the graph.
    * @post Nodes touching edge are erased as neighbors.
    * @post Number of active edges is reduced by 1. 
    * Complexity : no more than O(num_edges() + num_nodes()) */
  edge_iterator remove_edge(edge_iterator e_it) {
      remove_edge(*e_it);
      return this->edge_begin(); 
  }
 

  /** @brief Remove node @a n from the graph. 
    * @pre Node @a n is a valid node of the graph.
    * @return 1 to signal successful removal
    * @post Node @a n is removed from the graph. 
    * @post All edges incident to @a n are removed from the graph.
    * @post Number of active nodes is reduced by 1.
    * @post Adjacency vector of node is emptied.
    * Complexity : no more than O(num_nodes()) */
  size_type remove_node(const Node& n){
      unsigned oid = n.index();
      unsigned uid = ni2u_[oid];
      //Remove all incident edges
      unsigned i=0;
      while(!nodes_vec[uid].intadjnodes_.empty()){
          unsigned b_uid = nodes_vec[uid].intadjnodes_[i].adjnodeidx_;
          unsigned e_uid = nodes_vec[uid].intadjnodes_[i].adjeidx_;
          remove_edge(Edge(this, e_uid, uid, b_uid));
      }      
      //Swap and pop
      std::iter_swap(ni2u_.begin()+oid, ni2u_.end()-1);
      ni2u_.pop_back();
      //Set the replacement node's oid
      nodes_vec[ni2u_[oid]].onidx_ = oid;
      size_ = size_ - 1;
      return 1;
  }

 /** Remove node pointed to by @a n_it from the graph.
   * @brief Dereferences and calls remove_node on resulting node.
   * @pre Node pointed at is a valid node of the graph.
   * @return an iterator pointing to the first node.
   * @post Node is removed from the graph.
   * @post All edges incident are removed from the graph.
   * @post Adjacency vector of node is emptied.
   * @post Number of active nodes is reduced by 1. 
   * Complexity : no more than O(num_nodes()) */
  node_iterator remove_node(node_iterator n_it) {
      remove_node(*n_it);
      return this->node_begin(); 
  }

 




 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
