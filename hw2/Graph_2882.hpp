#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>
#include <cmath>
#include <unordered_set>
#include <iterator>
#include <utility>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V=int, typename E=int> // adding default values to templates
class Graph {
  // private:
  struct internal_node; // predefining structs
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node values. */
  using node_value_type = V;

  /** Type of edge values. */
  using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
    : gnodes_(), guid_(), next_uid_(0), gedges_(), 
    geuid_(), edge_check_(), next_euid_(0) {
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
  class Node: private totally_ordered<Node>  {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos;
    }

    /** Return this node's position which is modifiable. */
    Point& position(){
      return fetch().pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return grap_->guid_[fetch().uid];
    }

    /** Return this node's value. */
    node_value_type& value(){
      return fetch().value;
    }

    /** Return this node's const value. */
    const node_value_type& value() const{
      return fetch().value;
    }

    // Return the number of incident edges. */
    size_type degree() const{
      return fetch().connections.size();
    }

    // Start the incident iterator. */
    incident_iterator edge_begin () {
      return IncidentIterator(fetch().connections.begin(), *this, grap_);
    }

    // End the incident iterator. */
    incident_iterator edge_end() {
      return IncidentIterator(fetch().connections.end(), *this, grap_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      bool answer = (fetch().uid == n.uid_) and (this->grap_ == n.grap_);
      return answer;
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
      bool answer = (fetch().uid< n.uid_);
      return answer;
    }

   private:

    graph_type* grap_; // pointer back to the graph
    size_type uid_; // uid for every node
    
    Node(const graph_type* grap, size_type uid) // constructor
      : grap_(const_cast<graph_type*>(grap)), uid_(uid){
    }

    internal_node& fetch() const{ // function to get internal node from node
      return grap_->gnodes_[grap_->guid_.at(uid_)]; // hash map for constant time
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return gnodes_.size();
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

  Node add_node(const Point& position, const node_value_type& val=node_value_type()){
    internal_node new_node; // create new internal node
    new_node.pos = position;
    new_node.uid = next_uid_;
    new_node.value = val;
    gnodes_.push_back(new_node); // add it to graph containers
    guid_[next_uid_] = gnodes_.size() - 1;
    ++next_uid_;
    return Node(this, next_uid_-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return !(guid_.find(n.uid_) == guid_.end());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // assert(guid_.find(i) != guid_.end()); 
    size_type node_id = gnodes_[i].uid;
    return Node(this, node_id);        // Invalid node
  }

  /**
   * @brief Remove a given node from the graph.
   * 
   * @param a Node to be removed
   * @pre @a a is a node which can or cannot be in the graph
   * @post If old has_node(@a a), new num_nodes() == old num_nodes() - 1
   *       Else,                  new num_nodes() == old num_nodes()
   * @return size_type 1 if the node is found and removed, 0 otherwise
   * 
   * Complexity: O(a.degree())
   */
  size_type remove_node(const Node& a){
    Node& rem = const_cast<Node&>(a); // remove const to call incident iterator
    if(!has_node(a)){ // check if node exists
      return 0;
    }
    auto iit = rem.edge_begin(); // incident iterator
    for(; iit!= rem.edge_end(); ){
      Edge e_rem = *iit;
      remove_edge(e_rem); // remove all edges incident to node
      iit = rem.edge_begin(); // set it to begin since first element is erased 
    }

    size_type loc_erase = guid_[a.uid_]; // change index of last element
    guid_[gnodes_.back().uid] = loc_erase;
    std::iter_swap(gnodes_.begin()+loc_erase, gnodes_.end()-1); // swap and pop
    gnodes_.pop_back();

    return 1; // return 1 if node existed
  }

  /**
   * @brief Remove a given node from the graph.
   * 
   * @param n_it Node iterator pointing to the node to be removed
   * @pre @a n_it is a valid node iterator in the graph
   * @post If old has_node(@a *n_it), new num_nodes() == old num_nodes() - 1
   *       Else,                      new num_nodes() == old num_nodes()
   * @return node_iterator pointing to a node in the graph
   * 
   * Complexity: O((*n_it).degree())
   */
  node_iterator remove_node(node_iterator n_it){ // remove with iterator
    Node a = *n_it;
    ++n_it;
    remove_node(a); // call remove node using node
    return node_begin(); // returning an iterator that can be removed
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
  class Edge: private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return efetch().n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return efetch().n2;     // Invalid Node
    }

    /** Return the length of the edge. */
    double length() const{
      return efetch().elength;
    }

    /** Return this edge's value. */
    edge_value_type& value(){
      return efetch().value;
    }

    /** Return this edge's const value. */
    const edge_value_type& value() const{
      return efetch().value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.euid_==this->euid_) and (egrap_==e.egrap_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(egrap_==e.egrap_){
        return (this->euid_ < e.euid_); 
      }
      else{
        return egrap_<e.egrap_; // according to the test, edges from different graphs compare true vacuously
      }
    }

   private:
    
    graph_type* egrap_; // pointer back to graph
    size_type euid_; // eid for every edge
    
    Edge(const graph_type* egrap, size_type euid)
      : egrap_(const_cast<graph_type*>(egrap)), euid_(euid){
    }

    internal_edge& efetch() const{ // function to get internal edge from edge
      return egrap_->gedges_[egrap_->geuid_.at(euid_)]; // hash map for constant time
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return gedges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type edge_id = gedges_[i].eid;
    return Edge(this, edge_id);         // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::pair<size_type, size_type> check1 = {a.uid_, b.uid_}; // nodes can be in any order
    std::pair<size_type, size_type> check2 = {b.uid_, a.uid_};
    bool answer1 = (edge_check_.find(check1) != edge_check_.end());
    bool answer2 = (edge_check_.find(check2) != edge_check_.end());
    return (answer1 or answer2);
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
    std::pair<size_type, size_type> check1 = {a.uid_, b.uid_};
    std::pair<size_type, size_type> check2 = {b.uid_, a.uid_};

    if (edge_check_.find(check1) != edge_check_.end()){ // check if edge exists
      size_type edge_index = this->edge_check_[check1];
      return Edge(this, edge_index);
    }
    if (edge_check_.find(check2) != edge_check_.end()){
      size_type edge_index = this->edge_check_[check2];
      internal_edge new_edge; // if edge exists as (b, a)...
      new_edge.n1 = a; // ... we replace that edge by an edge of the form (a, b)...
      new_edge.n2 = b; // ... and return that edge
      new_edge.elength = norm(a.position()-b.position());
      new_edge.eid = edge_index;
      edge_check_.erase(check2); // erase edge (b, a) from map
      edge_check_[check1] = edge_index; // add edge (a,b) to the map
      gedges_[geuid_[edge_index]] = new_edge;
      return Edge(this, edge_index);
    }

    internal_edge new_edge; // if edge doesn't exist, add edge
    new_edge.n1 = a;
    gnodes_[guid_.at(a.uid_)].connections.push_back(b.uid_);
    new_edge.n2 = b;
    gnodes_[guid_.at(b.uid_)].connections.push_back(a.uid_);
    new_edge.eid = next_euid_;
    new_edge.elength = norm(a.position()-b.position());
    gedges_.push_back(new_edge);
    geuid_[next_euid_] = gedges_.size() - 1;
    edge_check_[check1] = next_euid_;
    ++next_euid_;
    return Edge(this, next_euid_-1); 

  }

  /**
   * @brief Remove a given edge from the graph using two nodes
   * 
   * @param a First node of the edge to be removed
   * @param b Second node of the edge to be removed
   * @pre @a a and @a b are valid nodes in the graph
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1
   *       Else,                        new num_edges() == old num_edges()
   * @return size_type 1 if the edge is found and removed, 0 otherwise
   * 
   * Complexity: O(log(num_edges())) due to finding in map
   */
  size_type remove_edge(const Node& a, const Node& b){

    if (!has_edge(a, b)){ // check if edge exists
      return 0;
    }
    // this block of code removes b from the vector of connections of a
    // and vice versa
    a.fetch().connections.erase(std::remove(a.fetch().connections.begin(), 
    a.fetch().connections.end(), b.uid_), a.fetch().connections.end());
    b.fetch().connections.erase(std::remove(b.fetch().connections.begin(), 
    b.fetch().connections.end(), a.uid_), b.fetch().connections.end());

    std::pair<size_type, size_type> check1 = {a.uid_, b.uid_};
    std::pair<size_type, size_type> check2 = {b.uid_, a.uid_};

    size_type eid_erase = 0; // get eid of edge to be removed
    if (edge_check_.find(check1) != edge_check_.end()){
      eid_erase = edge_check_[check1];
    }
    else if(edge_check_.find(check2) != edge_check_.end()){
      eid_erase = edge_check_[check2];
    }

    size_type loc_erase = geuid_[eid_erase]; // change index of last element
    geuid_[gedges_.back().eid] = loc_erase;
    std::iter_swap(gedges_.begin()+ loc_erase, gedges_.end()-1); // swap and pop
    gedges_.pop_back();

    edge_check_.erase(check1); // erase node pair from nodes-to-eid map
    edge_check_.erase(check2);

    return 1; // return 1 if edge existed
  }

  /**
   * @brief Remove a given edge from the graph.
   * 
   * @param e Edge to be removed from the graph
   * @pre @a e may or may not belong to the graph
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1
   *       Else,                        new num_edges() == old num_edges()
   * @return size_type 1 if the edge is found and removed, 0 otherwise
   * 
   * Complexity: O(log(num_edges())) due to finding in map
   */
  size_type remove_edge(const Edge& e){ // remove edge using edge
    Node a = e.node1();
    Node b = e.node2();
    return remove_edge(a, b); // call remove edge using nodes
  }

  /**
   * @brief Remove a given edge from the graph using edge iterator.
   * 
   * @param e_it Edge iterator pointing to the edge to be removed
   * @pre @a e_it is a valid edge iterator in the graph
   * @post If old has_edge(@a *e_it), new num_edges() == old num_edges() - 1
   *       Else,                      new num_edges() == old num_edges()
   * @return edge_iterator pointing to an edge in the graph
   * 
   * Complexity: O(log(num_edges()))
   */
  edge_iterator remove_edge(edge_iterator e_it) // using iterator
  {
    Edge e_rem = *e_it;
    ++e_it;
    remove_edge(e_rem); // call remove edge using edge
    return edge_begin(); // returning an iterator that can be removed
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() { // clear all containers and trackers
    gnodes_.clear();
    guid_.clear();
    next_uid_ = 0;
    gedges_.clear();
    geuid_.clear();
    edge_check_.clear();
    next_euid_ = 0;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }
    
    // Dereference operator
    Node operator*() const {
      if (nit_ == grap_->gnodes_.end()){ // return invalid node on dereferencing end
        return Node();
      }
      internal_node int_node = *nit_;
      return Node(grap_, int_node.uid);
    }

    // Increment operator
    NodeIterator& operator++(){
      ++nit_;
      return *this;
    }
    
    // Defines equality between two iterators
    bool operator==(const NodeIterator& nit) const{
      return ((nit.grap_ == this->grap_) and (nit.nit_ == this->nit_));
    }

    // Inequality operator not required due to inheritance

   private:
    friend class Graph;
    typename std::vector<internal_node>::const_iterator nit_; // iterator over internal nodes
    graph_type* grap_; // pointer back to the graph

    // Constructor
    NodeIterator(typename std::vector<internal_node>::const_iterator nit, const graph_type* grap)
    : nit_(nit), grap_(const_cast<graph_type*>(grap)) {}
  };

  // return an iterator pointing at the start of the vector of nodes
  node_iterator node_begin() const{
    return NodeIterator(gnodes_.begin(), this);
  }

  // return an iterator pointing at the end of the vector of nodes
  node_iterator node_end() const{
    return NodeIterator(gnodes_.end(), this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // Dereference operator
    Edge operator*() const{
      if(iit_ == iit_node_.fetch().connections.end()){ // return invalid edge on dereferencing end
        return Edge();
      }

      Node a = iit_node_;
      Node b = Node(grap_, *iit_);

      //return grap_->add_edge(a, b);
      std::pair<size_type, size_type> check1 = {a.uid_, b.uid_};
      std::pair<size_type, size_type> check2 = {b.uid_, a.uid_};

      if (grap_->edge_check_.find(check1) != grap_->edge_check_.end()){ // check if edge exists
        size_type edge_index = grap_->edge_check_[check1];
        return Edge(grap_, edge_index);
      }
      if (grap_->edge_check_.find(check2) != grap_->edge_check_.end()){
        size_type edge_index = grap_->edge_check_[check2];
        double old_elength = grap_->gedges_[grap_->geuid_[edge_index]].elength;
        internal_edge new_edge; // if edge exists as (b, a)...
        new_edge.n1 = a; // ... we replace that edge by an edge of the form (a, b)...
        new_edge.n2 = b; // ... and return that edge
        new_edge.eid = edge_index;
        new_edge.elength = old_elength;
        grap_->edge_check_.erase(check2); // erase edge (b, a) from map
        grap_->edge_check_[check1] = edge_index; // add edge (a,b) to the map
        grap_->gedges_[grap_->geuid_[edge_index]] = new_edge;
        return Edge(grap_, edge_index);
      }
      return Edge();
    }

    // Increment operator
    IncidentIterator& operator++(){
      ++iit_;
      return *this;
    }

    // Defines equality between two iterators
    bool operator==(const IncidentIterator& iit) const{
      return (iit.grap_ == grap_) and (iit.iit_node_==iit_node_) and (iit.iit_==iit_);
    }

    // Inequality operator not required due to inheritance

   private:
    friend class Graph;
    typename std::vector<size_type>::const_iterator iit_; // iterator over node_ids connected to spawning node
    Node iit_node_; // spawning node
    graph_type* grap_; // pointer back to the graph

    // Constructor
    IncidentIterator(typename std::vector<size_type>::const_iterator iit, Node iit_node, const graph_type* grap)
    : iit_(iit), iit_node_(iit_node), grap_(const_cast<graph_type*>(grap)) {}
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
    }

    // Dereference operator
    Edge operator*() const{
      if (eit_ == grap_->gedges_.end()){ // return invalid edge on dereferencing end
        return Edge();
      }
      internal_edge int_edge = *eit_;
      return Edge(grap_, int_edge.eid);
    }

    // Increment operator
    EdgeIterator& operator++(){
      ++eit_;
      return *this;
    }

    // Defines equality between two operators
    bool operator==(const EdgeIterator& eit) const{
      return ((eit.grap_ == this->grap_) and (eit.eit_ == this->eit_));
    }

   private:
    friend class Graph;
    typename std::vector<internal_edge>::const_iterator eit_; // iterator over internal edges
    graph_type* grap_; // pointer back to the graph

    // Constructor
    EdgeIterator(typename std::vector<internal_edge>::const_iterator eit, const graph_type* grap)
    : eit_(eit), grap_(const_cast<graph_type*>(grap)) {}
  };

  // return an iterator pointing at the start of the vector of edges
  edge_iterator edge_begin() const{
    return EdgeIterator(gedges_.begin(), this);
  }

  // return an iterator pointing at the end of the vector of edges
  edge_iterator edge_end() const{
    return EdgeIterator(gedges_.end(), this);
  }
  

 private:

  // for nodes

  struct internal_node{ // internal node for proxy setting
    Point pos;
    size_type uid;
    node_value_type value;
    std::vector<size_type> connections;
  };

  std::vector<internal_node> gnodes_; // vector to store nodes
  std::unordered_map<size_type, size_type> guid_; // mapping of uid to index in vec
  size_type next_uid_; // uid tracker
  
  // for edges

  struct internal_edge{ // internal edge for proxy setting
    Node n1;
    Node n2;
    size_type eid;
    double elength;
    edge_value_type value;
  };

  std::vector<internal_edge> gedges_; // vector to store edges
  std::unordered_map<size_type, size_type> geuid_; // mapping for edges
  std::map<std::pair<size_type,size_type>, size_type> edge_check_; // map from nodes to eid
  size_type next_euid_; // eid tracker

  Graph(const graph_type&) = delete;
  graph_type& operator=(const graph_type&) = delete;
};

#endif // CME212_GRAPH_HPP
