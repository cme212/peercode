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

template <typename V> // HW1
class Graph {
 private:

  // HW0
  struct node_internal;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V; //HW1

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
  /**
   * HW0: Every Graph starts with 100-size vector for node internals,
   * 100-size vector of Edges,
   * counters set to zero
   * indices for nodes and edges set to zero;
   * HW1: TODO: a 100-size vector of Nodes
   * deleted graph_size and edge_count because they are tracked by nuid and euid
   */
    : nint_vec_(0), edge_vec_(0), node_vec_(0),    
      next_nuid_(0), next_euid_(0) {
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
      : ng_ptr(nullptr), nuid(0) {
      // HW0:
      // Consider an invalid node to be one that has 
      // graph_node_ (pointer) == nullptr
      // and/or, nuid_ == 0;
      // HW1: renamed
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      Graph* my_graph = this->ng_ptr;
      return my_graph->nint_vec_[nuid].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nuid; // HW1 changed name to nuid
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    // HW1 TODID
    // QUESTION!! Why is this one not void? (if it's only assigning?)
    node_value_type& value() {
      Graph* my_graph = this->ng_ptr;
      return my_graph->nint_vec_[nuid].ivalue_;
    }
    
    // TODIDHW1!! QUESTION:: WHAT IS THE DIFFERENCE??
    const node_value_type& value() const {
      Graph* my_graph = this->ng_ptr;
      return my_graph->nint_vec_[nuid].ivalue_;
    }

    // returns the number of edges incident to this node
    size_type degree() const {
     return (*ng_ptr).adj_nuid_vec_[nuid].size();
    }

    IncidentIterator edge_begin() const {
      if (ng_ptr->adj_nuid_vec_[nuid].empty()) {
        return IncidentIterator(); // if you have no adjacents, return an invalid iterator
      } else {
        return IncidentIterator(ng_ptr, nuid, 0); // otherwise, start at zero
      }
    }

    IncidentIterator edge_end() const {
      if (ng_ptr->adj_nuid_vec_[nuid].empty()) {
        return IncidentIterator(); // if you have no adjacents, return an invalid iterator
      } else {
        return IncidentIterator(ng_ptr, nuid, ng_ptr->adj_nuid_vec_[nuid].size()); // o.w. return the last edge
      }
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool gr_eq = (this->ng_ptr == n.ng_ptr); // graphs are the same
      bool nuid_eq = (this->nuid == n.nuid); // unique id is the same
      return gr_eq && nuid_eq; // if both true, Nodes are the same
    } // HW1 : clunky but fine

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
        return this->nuid < n.nuid;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // pointer back to its Graph container // Step 8
    Graph* ng_ptr; // HW1 : changed name
    // this Node's unique ID
    size_type nuid; // HW1: changed name

    /** Private constructor */ // Step 9
    Node(const Graph* graph, size_type number)
        : ng_ptr (const_cast<Graph*>(graph)), nuid(number) {
        }
    
    /** Helper method to return the appropriate element
     * ng_ptr-> points to the graph's vector of node_internals
     * returns the the node_internal
     * that has the same uid as this Node
    */ // Step 10
    node_internal& fetch() const {
      // first get the address/pointer
        node_internal this_int_node = ng_ptr->nint_vec_[nuid];
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
    // HW1: updated to just next nuid, which already stores info on num Nodes
    return next_nuid_;
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
  Node add_node(const Point& p, const node_value_type& v = node_value_type()) {
    // HW0 and HW1
    // construct your new node_internal, add it to the nint_vec_
    node_internal this_int_node;
    this_int_node.position = p;
    this_int_node.ivalue_ = v;
    nint_vec_.push_back(this_int_node);

    // construct the corresponding Node
    Node new_node(this, next_nuid_);
    node_vec_.push_back(new_node);
    // add an empty vector in the corresponding adj_vec slot 
    // this is where you will bump nuid's of Nodes that connect through edges to this one
    std::vector<size_type> temp_adj_list(0); 
    adj_nuid_vec_.push_back(temp_adj_list);
    ++next_nuid_; // prepare for the next addition
    return new_node; // returns the Node constructed
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.ng_ptr == this;
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
    assert((i < size()) && (i >= 0));
    return node_vec_[i];
    // HW1 changed to use the vector of Nodes
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
      : eg_ptr(nullptr) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0
      return n1; // HW1 : changed for new Edge build
      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0
      return n2; // HW1 : changed for new build
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: 
      // equal edges match nodes either pairwise or opposite-wise
      bool check1 = (n1 == e.n1 && n2 == e.n2);
      bool check2 = (n1 == e.n2 && n2 == e.n1);
      return check1 || check2;
      // HW1 : updated with new build
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
      size_type te_n1 = n1.index(); // nuid of this first node
      size_type te_n2 = n2.index(); // nuid of this second node
      size_type oe_n1 = e.n1.index(); // nuid of other first node
      size_type oe_n2 = e.n2.index(); // nuid of other second node

      if (te_n1 + te_n2 == oe_n1 + oe_n2) {
        return std::min(te_n2, te_n1) < std::min(oe_n2, oe_n1);
      } else {
        return (te_n1 + te_n2) < (oe_n1 + oe_n2);
      } 
      return false; // quiet compiler warning
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0:
    Graph* eg_ptr; // know which Graph you belong to
    // HW1 : updated with two Nodes
    Node n1; 
    Node n2;

    // private constructor
    Edge(Graph* graph): eg_ptr(graph) {
      }
    // overloaded constructor
    
    Edge(const graph_type* graph, Node n1, Node n2)
      : eg_ptr(const_cast<graph_type*>(graph)),
        n1(n1), n2(n2) {
      } // if it ain't broke!!
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return next_euid_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < next_euid_);
    return edge_vec_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // iterates over a's list of adjacent nuids
    // if any of them match b's nuid, you have an edge
    for (size_type i = 0; i < a.degree(); i++) {
      if (adj_nuid_vec_[a.index()][i] == b.index()) { 
        return true;
      }
    }
    return false; // if no match, no edge
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

    // if an edge with those nodes exists, return it in the desired order
    if (has_edge(a, b)) {
      return Edge(this, a, b);
    }

    // if no edge exists, build one, add to vector
    Edge new_edge(this, a, b);
    edge_vec_.push_back(new_edge);
    next_euid_++;

    //update the adjacency lists
    adj_nuid_vec_[a.index()].push_back(b.index());
    adj_nuid_vec_[b.index()].push_back(a.index());

    // satisfy the function call
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    
    // wipe adjacency lists and Edges, reset count
    adj_nuid_vec_.clear();
    edge_vec_.clear();
    next_euid_ = 0;
    // wipe Nodes and reset count
    node_vec_ .clear();
    next_nuid_ = 0;
    // wipe internal nodes
    nint_vec_.clear();
    
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      node_itr_ptr = nullptr;
      node_itr_index = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    // increments the NodeIterator and returns the next NodeIterator
    // if such an action is valid
    NodeIterator operator++() {
      node_itr_index++;
      return (*this);
      /**
      if (node_itr_index < total_nodes) {
        node_itr_index++;
        return (*this);
      } else {
        return NodeIterator(ig_ptr, total_nodes);
      }*/
    }

    // returns true if two NodeIterator objects are pointing to the same index
    // of the same vector of nodes
    bool operator==(const NodeIterator& n_itr) const {
        return (node_itr_ptr == n_itr.node_itr_ptr) && (node_itr_index == n_itr.node_itr_index);
    }

    /**
    bool operator!=(const NodeIterator& n_itr) {
        return (this->node_iter_ptr != n_itr.node_itr_ptr) || (this->node_itr_index != n_itr.node_itr_index);
    }*/ // trying this out
    
    Node operator*() const { 
      assert(node_itr_index < ig_ptr->size()); // QUESTION!! Had strict equality, trying \leq
      return (*node_itr_ptr)[node_itr_index];
    }

    size_type my_nuid() {
      return node_itr_index;
    }

   private:
    friend class Graph;
    friend class EdgeIterator;
    // HW1 #2: YOUR CODE HERE
    Graph* ig_ptr;
    std::vector<Node>* node_itr_ptr;
    size_type node_itr_index;
    // size_type total_nodes; // having a set variable makes your NIter less flexible

    // private constructor that can be accessed by the Graph class
    // inputs a pointer to node vector and start index, given the number of nodes by the Graph variable next_nuid
    NodeIterator(const graph_type* graph, size_type index_start)
      : ig_ptr(const_cast<graph_type*>(graph)) {
      node_itr_ptr = &(ig_ptr->node_vec_);

      node_itr_index = std::min(ig_ptr->size(), index_start);
      /**CHANGE
      assert(index_start <= total_nodes); // had <= , made strict QUESTION!!
      node_itr_index = index_start;
      */
    }

  }; // END OF NODE ITERATOR CLASS


  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  NodeIterator node_end() const {
    return NodeIterator(this, next_nuid_); 
    // note that next_nuid is one past the last constructed node
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
    // An IncIter is invalid if it does not point to a Graph
    IncidentIterator() {
      inc_ptr = nullptr;
      nuid = 0;
      inc_itr_index = 0;
    }


    IncidentIterator& operator++() {
      inc_itr_index++;
      return (*this);
    }

    // Two IncIters are equal if they point to the same adj_vec and at the same index
    bool operator==(const IncidentIterator& some_itr) const {
      return (inc_ptr == some_itr.inc_ptr) && (inc_itr_index == some_itr.inc_itr_index)
        && (nuid == some_itr.nuid);
    }

    Edge operator*() const {
      // assert(inc_itr_index < num_conn_nodes);
      // size_type other_nuid = (*my_node_adj)[inc_itr_index];
      Node n1 = inc_ptr->node_vec_[nuid]; // retrieve the Node associated with this iterator
      unsigned other_nuid = inc_ptr->adj_nuid_vec_[nuid][inc_itr_index];
      Node n2 = inc_ptr->node_vec_[other_nuid]; // retrieve the Node at the other end of the existing Edge
      return Edge(inc_ptr, n1 ,n2); // return the Edge between them
    }


   private:
    friend class Graph;
    friend class EdgeIterator;
    // HW1 #3: YOUR CODE HERE
    graph_type* inc_ptr; // points to the Graph in question
    size_type nuid; // the nuid of the node which you are referencing from
    size_type inc_itr_index; // index into the vector of adjacent nuids

    // std::vector<size_type>* my_node_adj; // pointer to the vector of adjacent nuids
    // std::vector<Node>* nvec_ptr; // pointer to this Graph's vector of Nodes
    // size_type o_nuid; // the nuid of the node at the other end of the current Edge
    
    // size_type num_conn_nodes; // size of vector of adjacent nuids

    // private constructor accessible by Graph
    IncidentIterator(graph_type* graph, size_type uid, size_type index)
      : inc_ptr(graph), nuid(uid), inc_itr_index(index) {
        /**
        my_node_adj = &(inc_ptr->adj_nuid_vec_[nuid]);
        nvec_ptr = &(inc_ptr->node_vec_);
        num_conn_nodes = (*my_node_adj).size();
        if (inc_itr_index < num_conn_nodes) {
          o_nuid = (*my_node_adj)[inc_itr_index];
        }
        */
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
      eig_ptr = nullptr;
      currit_node = NodeIterator();
      endit_node = NodeIterator();
      currit_inc = IncidentIterator();
      endit_inc = IncidentIterator();
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    EdgeIterator& operator++() {
      
      // incident iter always increments
      ++currit_inc;
      if (currit_inc == endit_inc) { // if it's at the end, move to next Node
        ++currit_node;
        if (currit_node != endit_node) { // if your Node is valid, update incident inc
          currit_inc = (*currit_node).edge_begin();
          endit_inc = (*currit_node).edge_end();
        } else { // otherwise, invalidate it
          currit_inc = IncidentIterator();
          endit_inc = IncidentIterator();
        }
        return (*this);
      }
     return *this;
    }

    // two EdgeIterators are equal iff they point at the same Graph, 
    // and their NodeIterators and IncidentIterators index the same nodes
    bool operator==(const EdgeIterator& edgit) const {
      return (currit_node == edgit.currit_node) && (currit_inc == edgit.currit_inc);
    }

    Edge operator*() {
      return (*currit_inc);
    }



   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* eig_ptr;
    NodeIterator currit_node;
    NodeIterator endit_node;
    IncidentIterator currit_inc;
    IncidentIterator endit_inc;
    

    // private constructor
    // given a graph, node_index and inc_index,
    // constructs a NodeIterator at the Node from node_index
    // and an IncidentIterator at that same node indexed at inc_index in its adjacency list
    
    EdgeIterator(const graph_type* graph, size_type nuid) 
      : eig_ptr(const_cast<graph_type*>(graph)) {
        currit_node = NodeIterator(eig_ptr, nuid); // start where you're told
        endit_node = eig_ptr->node_end(); // always end at the end
        if (nuid < eig_ptr->size()){
          currit_inc = (*currit_node).edge_begin(); // first incident Edge of this Node
          endit_inc = (*currit_node).edge_end(); // last incident Edge of this Node
        }
        
      }
    
    
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  EdgeIterator edge_end() const {
    EdgeIterator end_edge = EdgeIterator(this, size()); // constructed at final Node
    end_edge.currit_inc = end_edge.endit_inc; // with its IncIter maxed out
    return end_edge;
  }

 private:
    friend class NodeIterator;
    friend class IncidentIterator;
    friend class EdgeIterator;
    
    struct node_internal { // Step 3
      public:
      Point position;
      node_value_type ivalue_;
    };
  

    std::vector<node_internal> nint_vec_; // the internal node objects 
    std::vector<std::vector<size_type> > adj_nuid_vec_; // for each Node, a list
    // of nuids of Nodes which it has and Edge with

    std::vector<Edge> edge_vec_; // vector of Edges 
    std::vector<Node> node_vec_; // vector of Nodes
    
    size_type next_nuid_; // counter for incrementing nodes
    size_type next_euid_; // counter for incrementing edges

    // disable copy and assignment of a Graph // Step 5
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
  
}; // end of Graph Class

#endif // CME212_GRAPH_HPP