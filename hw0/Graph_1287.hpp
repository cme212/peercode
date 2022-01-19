#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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
         unsigned intnidx_; //Index of node
         std::vector<node_adjacent> intadjnodes_; //Vector of adjacent nodes
         internal_node(Point intpoint, unsigned intnidx)
               : intpoint_(intpoint), intnidx_(intnidx){}
      };
     unsigned size_; //Number of nodes
     // Internal type of edges:
     struct internal_edge{
         unsigned inteidx_; //Index of edge
         unsigned intnode1idx_; //Index of node1
         unsigned intnode2idx_; //Index of node2
         internal_edge(unsigned inteidx, unsigned intnode1idx, unsigned intnode2idx)
             : inteidx_(inteidx), intnode1idx_(intnode1idx), intnode2idx_(intnode2idx){}
     };
     unsigned num_edges_; //Number of edges
     std::vector<internal_edge> edges_vec; //Vector of internal edge objects
     std::vector<internal_node> nodes_vec; //Vector of internal node objects
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
  class Node {
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

    /** Return this node's index, a number in the range [0, graph_size). */
    unsigned index() const {
      // HW0: CODE ADDED HERE
      // Return index attribute
        return nidx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

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
        return (this-> nidx_ < n.nidx_); 
        (void) n;           // Quiet compiler warning
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
     return size_;
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
  Node add_node(const Point& position) {
    // HW0: CODE ADDED HERE
    // Add a new instance of internal_node to the back of the vector
    nodes_vec.emplace_back(internal_node(position, size_)); 
    size_ = size_ + 1; // Increment number of nodes
    (void) position;      // Quiet compiler warning
    return Node(this, size_-1);        // Return node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: CODE ADDED HERE
    // Check if node is present in the vector of nodes
    unsigned n_idx = n.index();
    for(unsigned i=0; i<nodes_vec.size();i++){
        if(nodes_vec[i].intnidx_==n_idx){
            return true;
        }
    }
    (void) n;            // Quiet compiler warning
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
    (void) i;             // Quiet compiler warning
    return Node(this, i);        // Return node at index i
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
  class Edge {
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
      return (this-> eidx_ < e.eidx_);
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
        : egraph_(const_cast<Graph*>(graph)), eidx_(eidx), n1idx_(n1idx), n2idx_(n2idx){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: CODE ADDED HERE
      return num_edges_; //Return number of edges
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: CODE ADDED HERE
    (void) i;             // Quiet compiler warning
    unsigned node1idx = edges_vec[i].intnode1idx_;
    unsigned node2idx = edges_vec[i].intnode2idx_;
    return Edge(this, i, node1idx, node2idx); // Return Edge at i
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: CODE ADDED HERE
    unsigned a_idx = a.index();
    unsigned b_idx = b.index();
    // For each node in the corresponding adjacent nodes vector,
      // check to see if the index matches node b
    for(unsigned i=0; i< nodes_vec[a_idx].intadjnodes_.size();i++){
        if(nodes_vec[a_idx].intadjnodes_[i].adjnodeidx_ == b_idx){
            return true;   
    	} 
    }  
    (void) a; (void) b;   // Quiet compiler warning
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
    // Check if edge already exists
    unsigned a_idx = a.index();
    unsigned b_idx = b.index();
    // We loop again (instead of calling has_edge) to get the index
    for(unsigned i=0; i< nodes_vec[a_idx].intadjnodes_.size();i++){
        if(nodes_vec[a_idx].intadjnodes_[i].adjnodeidx_ == b_idx){
            return Edge(this, nodes_vec[a_idx].intadjnodes_[i].adjeidx_, a_idx, b_idx);  
        }
    }     
    // If the edge doesn't exist, add a new edge
    internal_node inode_a = nodes_vec[a_idx];
    internal_node inode_b = nodes_vec[b_idx];
    // Add two new instances of adjacent nodes (one for a, one for b)
    // Add these instances to the corresponding vector
    nodes_vec[a_idx].intadjnodes_.push_back(node_adjacent(num_edges_,b_idx));
    nodes_vec[b_idx].intadjnodes_.push_back(node_adjacent(num_edges_,a_idx));
    // Add one new instance of internal edge and add to vector
    edges_vec.push_back(internal_edge(num_edges_,a_idx,b_idx));
    // Increment number of edges
    num_edges_ = num_edges_ + 1;
    (void) a, (void) b;   // Quiet compiler warning
    return Edge(this, num_edges_-1,a_idx,b_idx); // Return New Edge
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: CODE ADDED HERE
     edges_vec.clear();
     nodes_vec.clear();
     size_ = 0;
     num_edges_ = 0;
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
