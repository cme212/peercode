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

template <typename V, typename E> // HW1
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
  using edge_value_type = E;

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
  Graph() 
  /**
   * 
   */
    : nint_vec_(0), node_vec_(0), 
      eint_vec_(0), edge_vec_(0), adj_euid_vec_(0),   
      next_nuid_(0), next_euid_(0) {
        //adj_uid_vec_(0), 
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
      : ng_ptr(nullptr), nuid_(0) {
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
      return my_graph->nint_vec_[nuid_].position;
    }

    /**
     * @brief public function that allows modification of a Node's position
     * @param[out] position non-constant referent to this Node's 
     *    member variable of type Point
     *  HW2 COMPLETE
     */
    Point& position() {
      return ng_ptr->nint_vec_[nuid_].position;
    }

    /**
     * @brief read-only public function that returns a Node's INITIAL position
     * @param[out] p0_ non-constant referent to this Node's 
     *    member variable of type Point
     *  HW2 COMPLETE
     */
    const Point& p0() const {
      return ng_ptr->nint_vec_[nuid_].p0_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return ng_ptr->nint_vec_[nuid_].ndx; // HW2 changed name to nuid_
    }

    size_type nuid() const {
      return nuid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    // HW1 TODID
    node_value_type& value() {
      Graph* my_graph = this->ng_ptr;
      return my_graph->nint_vec_[nuid_].ivalue_;
    }
    
    // TODIDHW1!! QUESTION:: WHAT IS THE DIFFERENCE??
    const node_value_type& value() const {
      Graph* my_graph = this->ng_ptr;
      return my_graph->nint_vec_[nuid_].ivalue_;
    }

    // returns the number of edges incident to this node
    size_type degree() const {
     // return (*ng_ptr).adj_uid_vec_[nuid_].size();
     return (*ng_ptr).adj_euid_vec_[nuid_].size();
    }

    IncidentIterator edge_begin() const {
      if (ng_ptr->adj_euid_vec_[nuid_].empty()) { //changed adj_uid_vec_ nuid -> euid
        return IncidentIterator(); // if you have no adjacents, return an invalid iterator
      } else {
        return IncidentIterator(ng_ptr, nuid_, 0); // otherwise, start at zero
      }
    }

    IncidentIterator edge_end() const {
       //changed adj_uid_vec_ nuid -> euid
        //return IncidentIterator(ng_ptr, nuid_, ng_ptr->adj_euid_vec_[nuid_].size()); // o.w. return the last edge
      
      
      if (ng_ptr->adj_euid_vec_[nuid_].empty()) { //changed adj_uid_vec_ nuid -> euid
        return IncidentIterator(); // if you have no adjacents, return an invalid iterator
      } else {    //changed adj_uid_vec_ nuid -> euid
        return IncidentIterator(ng_ptr, nuid_, ng_ptr->adj_euid_vec_[nuid_].size()); // o.w. return the last edge
      }
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      bool gr_eq = (this->ng_ptr == n.ng_ptr); // graphs are the same
      bool nuid_eq = (this->nuid_ == n.nuid_); // unique id is the same
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
        return nuid_ < n.nuid();
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // pointer back to its Graph container // Step 8
    Graph* ng_ptr; // HW1 : changed name
    // this Node's unique ID
    size_type nuid_; // HW1: changed name

    /** Private constructor */ // Step 9
    Node(const Graph* graph, size_type number)
        : ng_ptr (const_cast<Graph*>(graph)), nuid_(number) {
        }
    
    /** Helper method to return the appropriate element
     * ng_ptr-> points to the graph's vector of node_internals
     * returns the the node_internal
     * that has the same uid as this Node
    */ // Step 10
    node_internal& fetch() const {
      // first get the address/pointer
        node_internal this_int_node = ng_ptr->nint_vec_[nuid_];
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
    // HW1: updated to just next nuid_, which already stores info on num Nodes
    // HW2: updated to return the size of node_vec_
    return node_vec_.size();
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
   */ // Step 12
  Node add_node(const Point& p, const node_value_type& v = node_value_type()) {
    // HW0 and HW1
    // construct your new node_internal, add it to the nint_vec_
    node_internal this_int_node;
    this_int_node.position = p;
    this_int_node.p0_ = p;
    this_int_node.ivalue_ = v;
    this_int_node.ndx = node_vec_.size();
    
    nint_vec_.push_back(this_int_node);

    // construct the corresponding Node
    Node new_node(this, next_nuid_);
    node_vec_.push_back(new_node);
    // add an empty vector in the corresponding adj_vec slot 
    // this is where you will bump nuid's of Nodes that connect through edges to this one
    std::vector<size_type> temp_nadj_list(0);
    std::vector<size_type> temp_eadj_list(0);
    // adj_uid_vec_.push_back(temp_nadj_list); ////changed adj_uid_vec_ nuid -> uid
    adj_euid_vec_.push_back(temp_eadj_list);
    ++next_nuid_; // prepare for the next addition
    return new_node; // returns the Node constructed
  }

  /**
   * @brief removes Node n from the parent graph
   * 
   * @param[in] n the Node to be removed
   * @pre has_node(n) == true
   * @post has_node(n) == false
   * @post for every Edge e initially incident to n, has_edge(e) == false
   *
   * Complexity: O( num_nodes ) amortized operations (assuming the graph is sparse). 
   * The largest container that this function iterates over is, 
   * over all incident edges e, the max(e.node1().degree(), e.node2().degree())
   * 
   * @return size_type true (1) if successful; 
   *                   false (0) if n is not in graph
   */

  size_type remove_node(const node_type& n) {

    if (has_node(n)) {
      // first, remove all edges incident to the node
      size_type removed_edge_counter(0);
      size_type tdeg = n.degree(); // must be fixed beforehand, or will keep changing
      for (size_type i = 0; i < tdeg; i++) {
      size_type teuid = adj_euid_vec_[n.nuid()][i]; // the uid of this edge
      size_type tedx = eint_vec_[teuid].edx; // the index of this edge
      edge_type e = edge_vec_[tedx]; // the actuall edge
      removed_edge_counter += remove_edge(e);
      
    }

      /**
      size_type removed_edge_counter(0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        std::cout<<"Removing edge " << (*it).euid() << std::endl;
        
        removed_edge_counter += remove_edge(*it);
        std::cout<<"Removed edges: " << removed_edge_counter << std::endl;
      }*/

      size_type tndx = n.index();
      node_vec_[tndx] = node_vec_.back(); // replace spot with last node in vector
      size_type o_nuid = node_vec_[tndx].nuid(); // get the replacement Node's uid
      nint_vec_[o_nuid].ndx = tndx; // update replacement node's index
      node_vec_.pop_back();

      return 1;
    } else {
      return 0; // cannot remove a Node that is not there
    }
  }
  /**
   * @brief Iterator version of a node remover. Removes the
   * node that the node_iterator type is currently pointing to
   * 
   * @post num_node = previous num_nodes - 1, unless the graph
   * had no nodes before the call
   * @post num_edges = previous num_edges - (num of edges incident
   * to *n_it before the function call)
   * 
   * @param n_it
   * @return node_iterator the same iterator pointing to the same
   * "slot" of the node vector, but now with a different node 
   * (or no node, if it was the last node)
   * 
   * Complexity: O( num_nodes ) amortized operations (assuming the graph is sparse). 
   * The largest container that this function iterates over is, 
   * over all incident edges e, the max(e.node1().degree(), e.node2().degree())
   */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it != node_end()){
      (*n_it) = node_vec_.back();
      node_vec_.pop_back();
    }
    return n_it;
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
    assert((i < size()) && (i >= 0));
    return node_vec_[i];
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
    Edge() 
      : eg_ptr(nullptr) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return eg_ptr->eint_vec_[euid_].node1; 
      
    }

    /** Return the other node of this Edge */
    Node node2() const {
  
      return eg_ptr->eint_vec_[euid_].node2; // HW2: changed for new edge int build
    }

    /* returns a modifiable reference to this edge's value. Compl: O(1)*/
    edge_value_type& value() {
      return eg_ptr->eint_vec_[euid_].evalue_;
    }

    const edge_value_type& value() const {
      return eg_ptr->eint_vec_[euid_].evalue_;
    }

    size_type index() const {
      // HW0: YOUR CODE HERE
      return eg_ptr->eint_vec_[euid_].edx; // HW2 changed name to nuid_
    }

    size_type euid() const {
      return euid_;
    }


    // returns the euclidean distance between the two nodes
    double length() {
      Point diff = node1().position() - node2().position();
      return norm(diff); // euclidean distance
    }

    // flips the node1 and node2 nodes, 
    // useful when required to return an edge in a given order
    //HW2 : added for new edge int build
    void swap() {
      Node temp_node = eg_ptr->eint_vec_[euid_].node1;
      eg_ptr->eint_vec_[euid_].node1 = eg_ptr->eint_vec_[euid_].node2;
      eg_ptr->eint_vec_[euid_].node2 = temp_node;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * 
     * Two edges are equal if they point to the same graph
     * and have the same euid_
     */
    bool operator==(const Edge& e) const {
      return ((eg_ptr == e.eg_ptr) && (euid_ == e.euid()));
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
             // Quiet compiler warning
      //HW2 : updated for new edge int build
      if (*this == e) {
        return false; // if they are equal, neither is < the other
      } else if (eg_ptr == e.eg_ptr) {
        return euid_ < e.euid();
      } else {
        return eg_ptr < e.eg_ptr;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0:
    Graph* eg_ptr; // know which Graph you belong to
    size_type euid_; // this never changes

    // private constructor
    Edge(graph_type* graph): eg_ptr(graph) {
      }

    // overloaded constructor
    
    Edge(const graph_type* graph, size_type k)
      : eg_ptr(const_cast<graph_type*>(graph)),
        euid_(k) {
      } // new implementation
    
  }; // end of Edge class



  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_vec_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < edge_vec_.size());
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
      size_type teuid = adj_euid_vec_[a.nuid()][i]; // the uid of this edge
      size_type tedx = eint_vec_[teuid].edx; // the index of this edge
      edge_type e = edge_vec_[tedx]; // the actual edge
      if ((e.node1() == b) || (e.node2() == b)) { // if either node matches b, the edge exits
        return true;
      }
    }

    return false; // if no match, no edge

    /* Might come back to this version later
    for (auto it = edge_begin(); it!= edge_end(); ++it) {
      bool check1 = (((*it).node1() == a) && ((*it).node2() == b));
      bool check2 = (((*it).node1() == b) && ((*it).node2() == a));
      if (check1 || check2) {
        return true;
      }
    }*/
    /**
    for (size_type i = 0; i < a.degree(); i++) {
      if (adj_uid_vec_[a.nuid()][i] == ) { // changed
        return true;
      }
    }*/
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
    
    // HW2 implementation::
    for (size_type i = 0; i < a.degree(); i++) {
      size_type teuid = adj_euid_vec_[a.nuid_][i]; // get the index of an incident edge
      size_type eid = eint_vec_[teuid].edx;
      bool check1 = (edge_vec_[eid].node1() == a && edge_vec_[eid].node2() == b);
      bool check2 = (edge_vec_[eid].node1() == b && edge_vec_[eid].node2() == a);
      if (check1 || check2) { // if you have a match, check for proper order
        if (edge_vec_[eid].node1() == a) {
          return edge_vec_[eid];
        } else { // if not in desired order, swap then return
          edge_vec_[eid].swap();
          return edge_vec_[eid];
        }
      } // so now you are returning the actual edge, not an expedient version of it
    }
    

    // if no edge exists, build and internal edge, populate
    edge_internal nEi;
    nEi.node1 = a;
    nEi.node2 = b;
 
    
    nEi.edx = edge_vec_.size(); // the next index is the size of the vector of active edges
    eint_vec_.push_back(nEi);

    Edge nEd(this, next_euid_); // this is always unique
    edge_vec_.push_back(nEd); // at index edx, the edge has euid next_euid_

    //update the adjacency lists

    // nodes get each other's nuid's // don't need this with new implementation

    // nodes both get the same euid update
    adj_euid_vec_[a.nuid()].push_back(next_euid_);
    adj_euid_vec_[b.nuid()].push_back(next_euid_);
    

    // satisfy the function call
    next_euid_++;
    return nEd;
  }

  /**
   * @brief removes Edge e from the parent graph
   * 
   * @param[in] e the Edge to be removed
   * @pre has_edge(e) == true
   * @post has_edge(e) == false
   * @post both Nodes of this Edge have been unassigned from the common Edge
   * 
   * @return size_type true (1) if successful; 
   *                   false (0) if e is not in graph
   * 
   *   * Complexity: O( max(e.node1().degree(), e.node2().degree()) ) amortized operations. 
   * The largest container that this function iterates over is
   * the larger of the two adj_euid_vec_
   */
  size_type remove_edge(const Edge& e) {
    node_type a = e.node1();
    node_type b = e.node2();
    if (!has_edge(a, b)) {
      return 0; // if this edge does not exist in the graph
    } else {
      // remove this Edge's euid from each of the Node's adjacency lists
      for (auto it = adj_euid_vec_[a.nuid()].begin(); it != adj_euid_vec_[a.nuid()].end(); it++) {
        if ((*it) == e.euid()) { 
          (*it) = adj_euid_vec_[a.nuid()].back(); // replace with the last euid
          adj_euid_vec_[a.nuid()].pop_back(); // remove the last euid (which is double now)
          break;
        }
      }

      for (auto it = adj_euid_vec_[b.nuid()].begin(); it != adj_euid_vec_[b.nuid()].end(); it++) {
        if ((*it) == e.euid()) { 
          (*it) = adj_euid_vec_[b.nuid()].back(); // replace with the last euid
          adj_euid_vec_[b.nuid()].pop_back(); // remove the last euid (which is double now)
          break;
        }
      }

      size_type tedx = e.index(); // get this Edge's index
      edge_vec_[tedx] = edge_vec_.back(); // replace current index spot with edge at end of vec
      size_type o_euid = edge_vec_[tedx].euid(); // get the euid of the replacement Edge
      eint_vec_[o_euid].edx = tedx; // update replacement internal_edge's index 
      edge_vec_.pop_back();
      return 1;
    }
  }

  /**
   * @brief Two-node remove edge: removes the edge between two given nodes, if 
   * such and edges exists. Updates the node adjacency lists to ensure that
   * they do not think they are still connected.
   * 
   * @param[in] n1
   * @param[in] n2 The two nodes at opposite ends of the edge to be removed 
   * @return size_type true (1) if successful; 
   *                   false (0) if e is not in graph
   * 
   * * Complexity: O( max(e.node1().degree(), e.node2().degree()) ) amortized operations. 
   * The largest container that this function iterates over is
   * the larger of the two adj_euid_vec_
   */
  size_type remove_edge(const node_type& n1, const node_type& n2) {
    //std::cout<<"Calling remove_edge of two nodes " <<std::endl;

      //std::cout<<"two-node remove edge searching" << std::endl;
      for (size_type i = 0; i < n1.degree(); i++) {
      size_type teuid = adj_euid_vec_[n1.nuid()][i]; // the uid of this edge
      size_type tedx = eint_vec_[teuid].edx; // the index of this edge
      edge_type e = edge_vec_[tedx]; // the actual edge
      if ((e.node1() == n2) || (e.node2() == n2)) {
        //std::cout<<"two-node removing edge " << e.euid() << std::endl;
        return remove_edge(e);
      } 
    }
    
    /*
    for (auto it = n1.edge_begin(); it != n1.edge_end(); ++it){
      if ((*it).node2() == n2) { // looking for the Edge whose n2 matches
        remove_edge(*it);
        std::cout<<"two-node remove_edge successful."<<std::endl;
        std::cout<<"Graph num edges is " << num_edges() << std::endl;
        return 1;  
      }
    }*/
    // std::cout<<"No edge to remove" <<std::endl;
    return 0; // if you can't find the Edge
  }

  /**
   * @brief Iterator version of an edge remover. Removes the
   * edge that the edge_iterator type is currently pointing to
   * 
   * @post num_edges = previous num_edges - 1, unless the graph
   * had no edges before the call
   * 
   * @param e_it
   * @return edge_iterator the same iterator pointing to the same
   * "slot" of the edge vector, but now with a different edge 
   * (or no edge, if it was the last edge)
   * 
   * Complexity: O( max(e.node1().degree(), e.node2().degree()) ) amortized operations. 
   * The largest container that this function iterates over is
   * the larger of the two adj_euid_vec_
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end()) {
      size_type oc = remove_edge(*e_it);
      //() = edge_vec_.back();
      //edge_vec_.pop_back();
    }
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    
    // wipe adjacency lists and Edges, reset count
    // adj_uid_vec_.clear(); //changed nuid -> uid // might come back to this
    adj_euid_vec_.clear();
    edge_vec_.clear();
    next_euid_ = 0;

    // wipe Nodes and reset count
    node_vec_ .clear();
    next_nuid_ = 0;
    // wipe internal nodes and edges
    nint_vec_.clear();
    eint_vec_.clear();
    
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

    // increments the NodeIterator and returns the next NodeIterator
    // validity is ensured by the caller of the iterator
    NodeIterator operator++() {
      node_itr_index++;
      return (*this);
    }

    // returns true if two NodeIterator objects are pointing to the same index
    // of the same vector of nodes
    bool operator==(const NodeIterator& n_itr) const {
        return (node_itr_ptr == n_itr.node_itr_ptr) && (node_itr_index == n_itr.node_itr_index);
    }
    
    Node operator*() const { 
      assert(node_itr_index < ig_ptr->size());
      return (*node_itr_ptr)[node_itr_index];
    }

    size_type my_nuid() const {
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
      
    }

  }; // END OF NODE ITERATOR CLASS


  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  NodeIterator node_end() const {
    return NodeIterator(this, size()); // CHANGED THISSSSS 1/31/22
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

      size_type teuid = inc_ptr->adj_euid_vec_[nuid][inc_itr_index]; // the euid of edge we're on
      size_type eind = inc_ptr->eint_vec_[teuid].edx; // the index of the edge that we're on
      size_type tndx = inc_ptr->nint_vec_[nuid].ndx; // the index of our center node
      
      Node n = inc_ptr->node_vec_[tndx]; // the center Node of this iterator
      if (n == inc_ptr->edge_vec_[eind].node1()) {
        return inc_ptr->edge_vec_[eind];
      } else {
        inc_ptr->edge_vec_[eind].swap();
        return inc_ptr->edge_vec_[eind];
      }
      
    } 


   private:
    friend class Graph;
    friend class EdgeIterator;
    // HW1 #3: YOUR CODE HERE
    graph_type* inc_ptr; // points to the Graph in question
    size_type nuid; // the nuid of the node which you are referencing from
    size_type inc_itr_index; // index into the vector of adjacent nuids

  

    // private constructor accessible by Graph
    IncidentIterator(graph_type* graph, size_type uid, size_type index)
      : inc_ptr(graph), nuid(uid), inc_itr_index(index) {
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
      //endit_node = NodeIterator();
      currit_inc = IncidentIterator();
      //endit_inc = IncidentIterator();
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    EdgeIterator& operator++() {
      
      // incident iter always increments
      ++currit_inc;
      if (currit_inc == (*currit_node).edge_end()) { // if it's at the end, move to next Node
        ++currit_node;
        if (currit_node != eig_ptr->node_end()) { // if your Node is valid, update incident inc
          currit_inc = (*currit_node).edge_begin();
        } else { // otherwise, invalidate it
          currit_inc = IncidentIterator();
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
      // std::cout<<"Edge * works " << std::endl;
      return (*currit_inc);
    }



   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* eig_ptr;
    NodeIterator currit_node;
    IncidentIterator currit_inc;
    

    // private constructor
    // given a graph, node_index and inc_index,
    // constructs a NodeIterator at the Node from node_index
    // and an IncidentIterator at that same node indexed at inc_index in its adjacency list
    
    EdgeIterator(const graph_type* graph, size_type idx) 
      : eig_ptr(const_cast<graph_type*>(graph)) {
        currit_node = NodeIterator(eig_ptr, idx); // start where you're told
        if (idx < eig_ptr->size()){
          currit_inc = (*currit_node).edge_begin(); // first incident Edge of this Node
        } 
      }
  };

/** // SAVED EDGEITERATOR
 * class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    // Construct an invalid EdgeIterator. 
    EdgeIterator() {
      eig_ptr = nullptr;
      currit_node = NodeIterator();
      endit_node = NodeIterator();
      currit_inc = IncidentIterator();
      endit_inc = IncidentIterator();
    }


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
      // std::cout<<"Edge * works " << std::endl;
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
    
    EdgeIterator(const graph_type* graph, size_type idx) 
      : eig_ptr(const_cast<graph_type*>(graph)) {
        currit_node = NodeIterator(eig_ptr, idx); // start where you're told
        endit_node = eig_ptr->node_end(); // always end at the end
        if (idx < eig_ptr->size()){
          currit_inc = (*currit_node).edge_begin(); // first incident Edge of this Node
          endit_inc = (*currit_node).edge_end(); // last incident Edge of this Node
        } 
      }
  }; */


  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  EdgeIterator edge_end() const {
    EdgeIterator end_edge = EdgeIterator(this, size()); // constructed at one past final Node
    end_edge.currit_inc = IncidentIterator(); // with its IncIter invalidated
    return end_edge;
  }

 private:
    friend class NodeIterator;
    friend class IncidentIterator;
    friend class EdgeIterator;
    
    struct node_internal {
      public:
      Point position;
      Point p0_;
      node_value_type ivalue_;
      size_type ndx;
    };

    struct edge_internal {
      public:
      Node node1;
      Node node2;
      size_type edx;
      edge_value_type evalue_;
    };
  

    std::vector<node_internal> nint_vec_; // the internal node objects, indexed by nuid_
    std::vector<Node> node_vec_; // vector of Nodes, indexed by ndx
    // if you want to revert, check for all instances of adj_uid_vec_
    //std::vector<std::vector<size_type> > adj_uid_vec_; // vector indexed by nuid_ of lists of adjacent nuid_'s

    std::vector<edge_internal> eint_vec_; // internal edge objects, indexed by euid
    std::vector<Edge> edge_vec_; // vector of Edges, indexed by edx
    std::vector<std::vector<size_type> > adj_euid_vec_;
    
    
    size_type next_nuid_; // counter for incrementing nodes
    size_type next_euid_; // counter for incrementing edges

    // disable copy and assignment of a Graph // Step 5
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
  
}; // end of Graph Class

#endif // CME212_GRAPH_HPP