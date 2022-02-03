#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

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
 private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  struct adj_node;
  struct internal_node;
  struct internal_edge;


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

  /** Type of node value. */
  using node_value_type = V;

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
  Graph() :
    node_num(0),
    node_uid_num(0),
    edge_num(0),
    edge_uid_num(0),
    nodes_container(),
    edges_container(),
    current_nodes(),
    current_edges()
    {} 

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
    Node() {}
    // a proxy node that contains the pointer to the graph and its proxy id
    Node(const graph_type* graph, size_type id) 
      : graph_(const_cast<graph_type*>(graph)), id_(id) {
        assert(id <= graph->node_num-1);
      }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_container[graph_->current_nodes[id_]].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    /** Return the value of this node */
    node_value_type& value() {
      return graph_->nodes_container.at(graph_->current_nodes.at(id_)).value;
    }
    const node_value_type& value() const {
      return graph_->nodes_container.at(graph_->current_nodes.at(id_)).value;
    }

    /** Return the number of incidents
     */
    size_type degree() const{
      size_type uid = graph_->current_nodes.at(id_);
      return graph_->nodes_container.at(uid).connecting_nodes.size();
    }

    /** default incident iterator */
    incident_iterator edge_begin() const{
      return incident_iterator(this,0);
    }

    incident_iterator edge_end() const{
      return incident_iterator(this,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(n.id_ == id_ && n.graph_ == graph_){
        return true; 
        // equal if same index and graph pointers point to the same address
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
      if(n.id_<id_){
        return true; // first evaluate on uid
      }
      if(n.id_==id_){
        if (n.graph_ < graph_){return true;} // second evaluate on graph pointer
      }
      return false; 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type id_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_num;
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
  Node add_node(const Point& position, 
  const node_value_type& val = node_value_type()) {
    // add internal node struct to container
    nodes_container[node_uid_num] = internal_node(position,node_num,val);
    // add node uid to the proxy list
    current_nodes.push_back(node_uid_num);
    ++node_uid_num;
    ++node_num;
    return Node(this, node_num-1);      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(n.graph_ == this && n.id_<node_num){
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
    assert((0<=i)&&(i<node_num));
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}
    // edge is also a proxy similar to node
    Edge(const graph_type* graph, size_type id) 
      : graph_(const_cast<graph_type*>(graph)), 
      id_(id) {}

    /** Return a node of this Edge */
    Node node1() const {
      size_type uid1 = 
        graph_->edges_container[graph_->current_edges[id_]].node1_uid;
      size_type pid1 = graph_->nodes_container[uid1].proxy_id;
      return Node(graph_, pid1);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type uid2 = 
        graph_->edges_container[graph_->current_edges[id_]].node2_uid;
      size_type pid2 = graph_->nodes_container[uid2].proxy_id;
      return Node(graph_, pid2);        
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(e.id_ == id_ && e.graph_ == graph_){
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
      if(e.id_<id_){
        return true; // first evaluate on uid
      }
      if(e.id_==id_){
        if (e.graph_ < graph_){return true;} // second evaluate on graph pointer
      }
      return false; 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type id_; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<edge_num);
    return Edge(this, i);    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // get node a's connecting nodes from node a's internal node struct
    std::vector<adj_node> node_a_link = 
    nodes_container.at(current_nodes.at(a.id_)).connecting_nodes;
    size_type b_uid = current_nodes.at(b.id_);
    // check whether b is a's connecting node
    for(size_type i=0; i<node_a_link.size(); i++){
      if(node_a_link[i].node2_uid == b_uid){
        return true;
      }
    }
    return false;
  }

  /** find the edge uid connecting the two nodes, helper function for add_edge
   * @pre @a a and @a b are valid nodes of this graph
   * @return -1 if edge does not exist, otherwise return the edge uid
   */
  int check_edge_uid(const Node& a, const Node& b) const {
    // get node a's connecting nodes from node a's internal node struct
    std::vector<adj_node> node_a_link = 
    nodes_container.at(current_nodes.at(a.id_)).connecting_nodes;
    size_type b_uid = current_nodes.at(b.id_);
    // check whether b is a's connecting node
    for(size_type i=0; i<node_a_link.size(); i++){
      if(node_a_link[i].node2_uid == b_uid){
        return (int) node_a_link[i].edge_uid;
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
  Edge add_edge(const Node& a, const Node& b) {
    size_type a_uid = current_nodes.at(a.id_);
    size_type b_uid = current_nodes.at(b.id_);

    // check if the edge already exist
    int check = check_edge_uid(a,b);
    if(check != -1){
      // find the proxy id of this existing edge
      size_type edge_pid = edges_container[(size_type) check].proxy_id;
      return Edge(this, edge_pid);
    }
    // add internal edge object to edge container
    edges_container[edge_uid_num] = internal_edge(edge_num, a_uid, b_uid);
    // add the uid of this new edge to the vector of current edges
    current_edges.push_back(edge_uid_num);
    // add connecting node to the internal nodes struct in nodes container
    nodes_container[a_uid].connecting_nodes.push_back(
      adj_node(b_uid, edge_uid_num));
    nodes_container[b_uid].connecting_nodes.push_back(
      adj_node(a_uid, edge_uid_num));

    ++edge_num;
    ++edge_uid_num;
  
    return Edge(this, edge_num-1);    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_container.clear();
    current_nodes.clear();
    edges_container.clear();
    current_edges.clear();
    node_num = 0;
    node_uid_num = 0;
    edge_num = 0;
    edge_uid_num = 0;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Construct a valid iterator */
    NodeIterator(const graph_type* graph, size_type id) : 
      graph_(const_cast<graph_type*>(graph)), 
      id_(id) {
        // make sure the constructed iterator is valid
        assert(id <= graph->node_num);  
      }

    bool operator==(const NodeIterator& iter2) const{
      if((graph_ == iter2.graph_)&&(id_ == iter2.id_)){
        return true;
      }
      return false;
    }

    Node operator*() const{
      return Node(graph_, id_);
    }

    NodeIterator& operator++(){
      id_++;

      return *this;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type id_;
  };

  node_iterator node_begin() const{
    return NodeIterator(this,0);
  };
  // node_iterator node_end() const
  node_iterator node_end() const{
    return NodeIterator(this,this->node_num);
  };


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
    IncidentIterator(const node_type* node, size_type id) :
      node_(const_cast<node_type*>(node)),
      id_(id) {
        // make sure the constructed iterator is valid
        assert(id <= node->degree());  
      }

    Edge operator*() const{
      size_type node_uid = node_->graph_->current_nodes.at(node_->id_);
      size_type edge_uid = 
        node_->graph_->nodes_container.at(node_uid).connecting_nodes[id_].edge_uid;
      size_type edge_pid = node_->graph_->edges_container.at(edge_uid).proxy_id;
      // if node1 of edge is not node_, swap node1 and node2 in the internal edge
      Edge current_e = Edge(node_->graph_, edge_pid);
      Node node1 = current_e.node1();
      if(node1 != *node_){
        node_->graph_->edges_container.at(edge_uid).node1_uid = 
          node_->graph_->current_nodes.at(current_e.node2().id_);
        node_->graph_->edges_container.at(edge_uid).node2_uid = 
          node_->graph_->current_nodes.at(node1.id_);
      }
      return Edge(node_->graph_, edge_pid);
    }

    IncidentIterator& operator++(){
      id_++;
      return *this;
    }

    bool operator==(const IncidentIterator& iter2) const{
      // equal if two iterators have the same id and the nodes are equal
      if((id_ == iter2.id_)&&(*node_ == *iter2.node_)){
        return true;
      }
      return false;

    }

   private:
    friend class Graph;
    node_type* node_; // pointer to proxy node
    size_type id_; // index in the node's connecting_nodes vector
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    /** Construct a valid iterator */
    EdgeIterator(const graph_type* graph, size_type id) : 
      graph_(const_cast<graph_type*>(graph)), 
      id_(id) {
        // make sure the constructed iterator is valid
        assert(id <= graph->edge_num);  
      }

    bool operator==(const EdgeIterator& iter2) const{
      if((graph_ == iter2.graph_)&&(id_ == iter2.id_)){
        return true;
      }
      return false;
    }

    Edge operator*() const{
      return Edge(graph_, id_);
    }

    EdgeIterator& operator++(){
      id_++;
      return *this;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type id_;
  };

  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  };
  // node_iterator node_end() const
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->edge_num);
  };

 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // internal type for node adjacency, used in internal_node
  struct adj_node{
    size_type node2_uid;
    size_type edge_uid;

    adj_node(){}
    adj_node(size_type node_id, size_type edge_id): 
      node2_uid(node_id), edge_uid(edge_id) {}
  };

  // internal type for nodes
  struct internal_node{
    Point position;
    // vector of adjacent nodes
    std::vector<adj_node> connecting_nodes;
    size_type proxy_id;
    node_value_type value;
    internal_node(){}
    internal_node(const Point& point, size_type id, const node_value_type& val): 
      position(point), proxy_id(id), value(val) {}
  };

  // internal type for edges
  struct internal_edge{
    size_type proxy_id;
    size_type node1_uid;
    size_type node2_uid;
    internal_edge(){}
    internal_edge(size_type pid, size_type id1, size_type id2) : 
      proxy_id(pid), node1_uid(id1), node2_uid(id2) {}
  };

  // variables to keep track of the valid & total number of nodes and edges
  size_type node_num;
  size_type node_uid_num;
  size_type edge_num;
  size_type edge_uid_num;

  // containers to store all nodes and edges
  // maps of uid : internal struct
  std::unordered_map<size_type,internal_node> nodes_container;
  std::unordered_map<size_type,internal_edge> edges_container;

  // containers to store the uid of all the currently valid nodes and edges
  // vector with length = node_num or edge_num
  // values are between [0, node_uid_num) or [0, edge_uid_num)
  std::vector<size_type> current_nodes;
  std::vector<size_type> current_edges;

};

#endif // CME212_GRAPH_HPP
