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
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct node_internal; // Step 1

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
  using size_type = unsigned; // Step 2 : set or verify

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph.
   * starting with a vector size 100 of nullptrs 
   * this will fill with pointers to node_internals
   * that correspond to the uid of a Node
   * where the uid matches the index of ptr in the vector 
   */
  Graph() // Step 6, add these
  // HW0: YOUR CODE HERE
  /**
   * Every Graph starts with 100-size vector for node internals,
   * 100-size vector of Edges,
   * counters set to zero
   * indices for nodes and edges set to zero;
   */
    : node_ints_(100), graph_size_(0), next_nuid_(0),
     edges_(100),  edge_count_(0), next_euid_(0) {
  }

  /** Default destructor */
  ~Graph() = default; // Step 7, TODO!!

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
      : graph_node_(nullptr), nuid_(0) {
      // HW0: YOUR CODE HERE
      // Consider an invalid node to be one that has 
      // graph_node_ (pointer) == nullptr
      // and/or, nuid_ == 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      Graph* my_graph = this->graph_node_;
      return my_graph->node_ints_[nuid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nuid_;
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
      // HW0: YOUR CODE HERE
      bool gr_eq = (this->graph_node_ == n.graph_node_); // graphs are the same
      bool nuid_eq = (this->nuid_ == n.nuid_); // unique id is the same
      return gr_eq && nuid_eq; // if both true, Nodes are the same
      // (void) n;          // Quiet compiler warning
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
      if (*this == n) {
        return false;
      } else {
        return this->nuid_ < n.nuid_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // pointer back to its Graph container // Step 8
    Graph* graph_node_; 
    // this Node's unique ID
    size_type nuid_;

    /** Private constructor */ // Step 9
    Node(const Graph* graph, size_type uid)
        : graph_node_ (const_cast<Graph*>(graph)), nuid_(uid) {

        }
    
    /** Helper method to return the appropriate element
     * graph_node_-> points to the vector of pointers to node_internals
     * returns the pointer corresponding to the node_internal
     * that has the same uid as the Node in question
    */ // Step 10
    node_internal& fetch() const {
      // first get the address/pointer
        node_internal this_int_node = graph_node_->node_ints_[nuid_];
        return this_int_node; // 
    }
    /**
     * This used a loop in the proxy example, but they didn't 
     * have the direct indexing
     */
  }; // End of Node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */ // Step 11
  size_type size() const {
    // HW0: YOUR CODE HERE
    return graph_size_;
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
   * 
   * Note: Assumes the Graph is constructed with an initial 
   * node_internal* vector of size 100
   */ // Step 12
  Node add_node(const Point& p) {
    // HW0: YOUR CODE HERE
    node_internal this_int_node;
    this_int_node.uid = next_nuid_;
    this_int_node.position = p;
    size_type vec_size = node_ints_.size();
    // check the filling status of node_ints_ vector:
    if (next_nuid_ > vec_size / 2) {
        // create a new vector, twice as big
        std::vector<node_internal> new_node_ints_(2 * vec_size);
        // copy // using graph_size because we only have that many elts
        for (size_type i = 0; i < graph_size_; i++) {
            new_node_ints_[i] = node_ints_[i];
        }
        // scrap the old vector, replace with the new, bigger vector
        // delete[] node_ints_; // TODO FIGURE THIS OUT!!
        node_ints_ = new_node_ints_;
    }
    // regardless, "add" the new Node by assigning
    node_ints_[next_nuid_] = this_int_node;
    ++graph_size_;
    ++next_nuid_; // prepare for the next addition
    return Node(this, next_nuid_ - 1); // returns the Node constructed
    // using the address of its parent Graph and the pointer to 
    // the node_internal corresponding to its uid
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_node_ == this;
    // return true if the graph_node_ pointer points to this Graph
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size());
    return Node(this, i);
    // have disabled copying so this is OK
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
      : g_edge_ptr_(nullptr), euid_(0),
      node_1(), node_2() {
      // HW0: YOUR CODE HERE
      // pointer to null and both Nodes invalid makes an invalid edge
      // TODO : Is this enough to " invalidate" an edge?
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return this->node_1;
      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return this->node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // equal edges match either pairwise or opposite-wise
      bool check1 = (this->node_1 == e.node_1 && this->node_2 == e.node_2);
      bool check2 = (this->node_1 == e.node_2 && this->node_2 == e.node_1);
      return check1 || check2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * ORDERING AS FOLLOWS:
     * Assume edges are not equal. 
     * Compare the pairwise sum of their node uids. 
     * If the sum is equal, compare their smallest node uids
     * 
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if (*this == e) return false; // if they are equal, neither is < the other
      size_type te_n1 = this->node_1.nuid_;
      size_type te_n2 = this->node_2.nuid_;
      size_type oe_n1 = e.node_1.nuid_;
      size_type oe_n2 = e.node_2.nuid_;

      if (te_n1 + te_n2 == oe_n1 + oe_n2) {
        return std::min(te_n2, te_n1) < std::min(oe_n2, oe_n1);
      } else {
        return (te_n1 + te_n2) < (oe_n1 + oe_n2);
      } 

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* g_edge_ptr_; // know which Graph you belong to
    size_type euid_; // this Edge's unique ID
    Node node_1; // TODO!! CAN I HAVE POINTERS TO NODES OR NEED NODES???
    Node node_2;

    // private constructor
    Edge(const Graph* graph, size_type euid)
      : g_edge_ptr_(const_cast<Graph*>(graph)), euid_(euid) {

      }

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_count_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return edges_[i];
    //return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type i = 0; i < num_edges(); i ++) {
        bool check1 = (edges_[i].node_1 == a && edges_[i].node_2 == b);
        bool check2 = (edges_[i].node_2 == a && edges_[i].node_1 == b);
        if (check1 || check2) return true;
    }
    return false;
    // POSSIBLY:: CONSTRUCT A TEMP NODE AND COMPARE TO ONES IN VECTOR? TODO!!
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
    // HW0: YOUR CODE HERE // TODO TODO!!!
    // if an edge with those nodes exists, return it
    for (size_type i = 0; i < num_edges(); i ++) {
        bool check1 = (edges_[i].node_1 == a && edges_[i].node_2 == b);
        bool check2 = (edges_[i].node_2 == a && edges_[i].node_1 == b);
        if (check1 || check2) return edges_[i];
    }
    // otherwise, construct a new edge
    Edge new_edge(this, next_euid_); // private constructor
    new_edge.node_1 = a;
    new_edge.node_2 = b;
    // check vector size before adding
    if (next_euid_ > edges_.size() / 2) {
      //copy existing edges into a vector twice as big
      std::vector<Edge> new_edges_(2 * edges_.size());
      for (size_type i = 0; i < num_edges(); i++) {
        new_edges_[i] = edges_[i];
      } // replace the old vector
    edges_ = new_edges_;
    }
    //regardless, add the new edge
    edges_[next_euid_] = new_edge;
    ++edge_count_;
    ++next_euid_;
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // replace the Edge vector with an empty one, reset counts
    edges_.clear();
    std::vector<Edge> new_edges_(100);
    edges_ = new_edges_;
    edge_count_ = 0;
    next_euid_ = 0;

    

    // replace node_internal vector with an empty one, reset counts
    node_ints_.clear();
    std::vector<node_internal> new_node_ints_(100);
    node_ints_ = new_node_ints_;
    graph_size_ = 0;
    next_nuid_ = 0;
    
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
    struct node_internal { // Step 3
        Point position;
        size_type uid;
    };
    
    // Step 4
    std::vector<node_internal> node_ints_; // the internal objects
    std::vector<Edge> edges_; // vector of Edges // TODO!! IS THIS LEGIT??
    size_type graph_size_; // the # of Nodes in this Graph
    size_type edge_count_; // the # of Edges in thie Graph
    size_type next_nuid_; // helper for incrementing nodes
    size_type next_euid_; // helper for incrementing edges

    // disable copy and assignment of a Graph // Step 5
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
  
}; // end of Graph Class

#endif // CME212_GRAPH_HPP