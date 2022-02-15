#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>
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
template <typename V = int, typename E = double>
class Graph {
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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  size_type graph_size_;
  size_type next_uid_;

  size_type num_edges_;
  size_type edge_next_uid_;
  struct node_proxy{
    Point point;
    size_type uid;
    node_value_type value;

    node_proxy();
    node_proxy(Point point, size_type uid) :
      point(point), uid(uid) {}
    node_proxy(Point point, size_type uid, node_value_type value) :
      point(point), uid(uid), value(value) {}
  };


  // hash map of node id: node_proxy
  std::unordered_map <size_type, node_proxy> nodes; 
  std::unordered_set<size_type> valid_nodes; // for checking existence
  std::vector<size_type> valid_nodes_list; // for iterating
  
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    graph_size_ = 0;
    next_uid_ = 0;
    num_edges_ = 0;
    edge_next_uid_ = 0;
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
    }

    /** Return this node's position. */
    Point& position() {
      return fetch().point;
    }

    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      size_type uid_to_find = fetch().uid;
      auto it = std::find(graph_->valid_nodes_list.begin(), graph_->valid_nodes_list.end(), uid_to_find);
      if(it != graph_->valid_nodes_list.end()){
        return std::distance(graph_->valid_nodes_list.begin(), it);
      }else{
        std::cout << "Error: Indexing into node that does not exist in graph" << std::endl;
        throw;
      }
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    // node_value_type& value();
    // functions for retrieving node values
    node_value_type& value() {
      return fetch().value;
    }
    // const node_value_type& value() const;

    const node_value_type& value() const{
      return fetch().value;
    }
    // size_type degree() const;
    size_type degree() const{
        size_type degree = graph_->adj_list.at(fetch().uid).size();
        return degree;
    }
    // incident_iterator edge_begin() const;
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, &(graph_->adj_list.at(fetch().uid)), this, 0);
    }

    // incident_iterator edge_end() const;

    incident_iterator edge_end() const{
      return IncidentIterator(graph_, &(graph_->adj_list.at(fetch().uid)), this, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      bool same_graph;
      bool same_index;

      same_graph = this->graph_ == n.graph_;
      same_index = this->uid_ == n.uid_;

      return same_graph && same_index;
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

      bool same_graph;
      bool same_index;

      same_graph = this->graph_ == n.graph_;
      same_index = this->uid_ == n.uid_;
      if(same_graph && same_index){
        return false;
      }else if(!same_graph){
        return &graph_ > &(n.graph_);
      }else if (same_graph && this->uid_ < n.uid_){
        return true;
      }else{
        return false;
      }
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type uid_;

    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    node_proxy& fetch() const{
      assert(graph_->valid_nodes.find(uid_) != graph_->valid_nodes.end());

      return graph_->nodes.at(uid_);

      assert(false);

    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
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
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning

    // create a new node
    node_proxy new_node(position, next_uid_, node_value);;
    nodes.insert(std::make_pair(next_uid_, new_node));
    valid_nodes.insert(next_uid_);
    valid_nodes_list.push_back(next_uid_);

    // add node as key to adj list
    std::vector<edge_proxy> empty_list;
    adj_list.insert(std::make_pair(next_uid_, empty_list));

    // increases graph size attr
    ++graph_size_;
    ++next_uid_;
    return Node(this, next_uid_-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return valid_nodes.find(n.uid_) != valid_nodes.end();
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
    size_type index = valid_nodes_list[i];
    return Node(this, index);        // Invalid node
  }
 private:

  struct edge_proxy{
    node_type node1;
    node_type node2;
    size_type uid;
    edge_value_type value;

    edge_proxy() {};
    edge_proxy(node_type node1, node_type node2, size_type uid) :
      node1(node1), node2(node2), uid(uid) {}
    edge_proxy(node_type node1, node_type node2, size_type uid, node_value_type value) :
      node1(node1), node2(node2), uid(uid), value(value) {}
  };
  std::unordered_map<size_type, edge_proxy> edges;
  std::unordered_set<size_type> valid_edges;
  std::vector<size_type> valid_edges_list;
  std::unordered_map<size_type, std::vector<edge_proxy> > adj_list;

 public:
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return fetch().node1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch().node2;      // Invalid Node
    }

    edge_value_type& value(){
      return fetch().value;
    }

    const edge_value_type& value() const{
      return fetch().value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      bool same_graph;
      bool same_edge;

      same_graph = this->graph_ == e.graph_;
      same_edge = this->uid_ == e.uid_;

      return same_graph && same_edge;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE

      bool same_graph;
      bool same_edge;

      same_graph = this->graph_ == e.graph_;
      same_edge = this->uid_ == e.uid_;

      if(same_graph && same_edge){
        return false;
      }else if(!same_graph){
        return &graph_ > &(e.graph_);
      }else if (same_graph && this->uid_ < e.uid_){
        return true;
      }else{
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type uid_;

    Edge(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid){

    }

    edge_proxy& fetch() const{
      assert(graph_->valid_edges.find(uid_) != graph_->valid_edges.end());
      return graph_->edges.at(uid_);
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    size_type index = valid_edges_list[i];
    return Edge(this, index);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    // if not same graph, then false
    if (a.graph_ != b.graph_){
      return false;
    }
    if (adj_list.find(a.uid_) != adj_list.end()){
      std::vector<edge_proxy> a_edges = adj_list.at(a.uid_);

      for(unsigned i = 0; i < a_edges.size(); i++){
        if((a_edges[i].node1 == a && a_edges[i].node2 == b) ||
          (a_edges[i].node1 == b && a_edges[i].node2 == a)){
            return true;
          }
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
    (void) a, (void) b;   // Quiet compiler warning

    bool new_edge = true;
    // check if this node already exists in graph
    if (has_edge(a, b)){
      new_edge = false;
    }

    // if new edge, input the new edge into map
    // else do nothing and return the edge already in the graph

    if(new_edge){
      // create a new node
      edge_proxy edge(a, b, edge_next_uid_);
      edges.insert(std::make_pair(edge_next_uid_, edge));
      valid_edges.insert(edge_next_uid_);
      valid_edges_list.push_back(edge_next_uid_);

      // add to adj list for both nodes
      adj_list.at(a.uid_).push_back(edge);
      adj_list.at(b.uid_).push_back(edge);

      // increases graph size attr
      ++num_edges_;
      ++edge_next_uid_;

      return Edge(this, edge_next_uid_-1);
    }else{
      size_type index = 0;
      std::vector<edge_proxy> a_edges = adj_list.at(a.uid_);

      for(auto e_prox_it = a_edges.begin(); e_prox_it != a_edges.end(); ++e_prox_it){
        if(((*e_prox_it).node1 == a && (*e_prox_it).node2 == b) ||
          ((*e_prox_it).node1 == b && (*e_prox_it).node2 == a)){
            index = (*e_prox_it).uid;
          }
      }

      edge_proxy edge(a,b, index);
      edges[index] = edge;
      return Edge(this, index);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    graph_size_ = 0;
    next_uid_ = 0;
    num_edges_ = 0;
    edge_next_uid_ = 0;

    nodes.clear();
    edges.clear();
    adj_list.clear();
    valid_nodes.clear();
    valid_nodes_list.clear();
    valid_edges.clear();
    valid_edges_list.clear();

    return;
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
    }

    NodeIterator (const graph_type* graph, size_type index) {
        this->graph = const_cast<graph_type*>(graph);
        this->index = index;
    }
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    value_type operator*() const{
      
      return graph->node(index);
    }
    // NodeIterator& operator++()
    NodeIterator operator++(){
      return NodeIterator(graph, ++index);
    }

    // bool operator==(const NodeIterator&) const
    bool operator==(const NodeIterator& iterator) const {
        return this->graph == iterator.graph && this->index == iterator.index;
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph;
    size_type index;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  }
  // node_iterator node_end() const

  node_iterator node_end() const{
      return NodeIterator(this, graph_size_);
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

    IncidentIterator(const graph_type* graph, std::vector<edge_proxy>* edges, const node_type* node, size_type index ) {
      this->graph = const_cast<graph_type*>(graph);
      this->edges = const_cast<std::vector<edge_proxy>* >(edges);
      this->node = const_cast<node_type*>(node);
      this->index = index;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      assert(index < (*edges).size());

      if((*edges)[index].node1 != *node){

        // swap the nodes in the edges *map* if wrong orientation
        node_type tmp = {static_cast<node_type&&>(graph->edges.at((*edges)[index].uid).node1)};
        graph->edges.at((*edges)[index].uid).node1 = static_cast<node_type&&>(graph->edges.at((*edges)[index].uid).node2);

        graph->edges.at((*edges)[index].uid).node2 = static_cast<node_type&&>(tmp);
      }
      return graph->edge((*edges)[index].uid);
    }
    // IncidentIterator& operator++()
    IncidentIterator operator++(){
      return IncidentIterator(graph, edges, node, ++index);
    }
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& inc_it) const{
      return this->edges == inc_it.edges && this->index == inc_it.index;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph;
    std::vector<edge_proxy>* edges;
    node_type* node;
    size_type index;
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

    EdgeIterator(const graph_type* graph, size_type index){
      this->graph = const_cast<graph_type*>(graph);
      this->index = index;
    }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      return graph->edge(index);
    }
    // EdgeIterator& operator++()
    EdgeIterator operator++(){
      return EdgeIterator(graph, ++index);
    }
    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& edge_it) const{
      return this->graph == edge_it.graph && this->index == edge_it.index;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph;
    size_type index;


  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }
  // edge_iterator edge_end() const
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges_);
  }

  /** remove an edge from the graph.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @pre has_edge( @a a, @a b) == true
   * @return the number of edges removed (0,1)
   * @post has_edge( @a a, @a b) == false
   * @post If old has_edge( @a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Complexity: No more than O(num_edges()) from linear search of valid edges
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(!has_edge(a,b)){
      std::cout << "ERROR: attempting to remove edge that doesn't exist" << std::endl;
      return 0;
    }

    // find right index to erase
    size_type index = 0;
    std::vector<edge_proxy> a_edges = adj_list.at(a.uid_);

    for(auto e_prox_it = a_edges.begin(); e_prox_it != a_edges.end(); ++e_prox_it){
      if(((*e_prox_it).node1 == a && (*e_prox_it).node2 == b) ||
        ((*e_prox_it).node1 == b && (*e_prox_it).node2 == a)){
          index = (*e_prox_it).uid;
        }
    }

    // erase in the valid edges set
    valid_edges.erase(index);

    // find in valid edges list
    auto it = std::find(valid_edges_list.begin(), valid_edges_list.end(), index);
    assert(it != valid_edges_list.end());

    // swap and pop from valid edges list
    *it = valid_edges_list.back();
    valid_edges_list.pop_back();

    // remove from adj lists
    std::vector<edge_proxy>* a_edges_pointer = &(adj_list.at(a.uid_));
    std::vector<edge_proxy>* b_edges_pointer = &(adj_list.at(b.uid_));
    for(auto edge_it = (*a_edges_pointer).begin(); edge_it != (*a_edges_pointer).end(); ){
      if((*edge_it).uid == index) {
        std::swap((*edge_it), ((*a_edges_pointer).back()));
        (*a_edges_pointer).pop_back();
        break;
      }else{
        ++edge_it;
      }
    }

    for(auto edge_it = (*b_edges_pointer).begin(); edge_it != (*b_edges_pointer).end(); ){
      if((*edge_it).uid == index) {
        std::swap((*edge_it), ((*b_edges_pointer).back()));
        (*b_edges_pointer).pop_back();
        break;
      }else{
        ++edge_it;
      }
    }

    // update num of edges
    --num_edges_;

    assert(valid_edges.size() == num_edges_);
    return 1;


  }

  /** remove an edge from the graph.
   * @pre Edge object @a e has valid nodes e.node1() and e.node2()
   * @pre has_edge(e.node1(), e.node2()) == true
   * @return the number of edges removed (0,1)
   * @post has_edge(e.node1(), e.node2()) == false
   * @post If old has_edge(e.node1(), e.node2()), new num_edges() == old num_edges() - 1.
   *       Else,                                  new num_edges() == old num_edges().
   *
   * Complexity: No more than O(num_edges()) from linear search of valid edges
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /** remove an edge from the graph.
   * @pre Edge iterator @a e_it can be dereferenced to a valid Edge object e
   * @pre Edge object e has valid nodes e.node1() and e.node2()
   * @pre has_edge(e.node1(), e.node2()) == true
   * @return an Edge iterator object that points to the beginning of valid edges
   * @post has_edge(e.node1(), e.node2()) == false
   * @post If old has_edge(e.node1(), e.node2()), new num_edges() == old num_edges() - 1.
   *       Else,                                  new num_edges() == old num_edges().
   *
   * Complexity: No more than O(num_edges()) from linear search of valid edges
   */
  edge_iterator remove_edge(edge_iterator e_it){
    size_type index = remove_edge((*e_it).node1(), (*e_it).node2());
    return this->edge_begin();
  }

  /** remove a node and all incident edges from the graph.
   * @pre @a n is a valid node of this graph
   * @pre has_node( @a n) == true
   * @return the number of nodes removed (0,1)
   * @post has_node( @a n) == false
   * @post If old has_node( @a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                   new num_nodes() == old num_nodes().
   * @post new n.degree() == 0.
   * @post new num_edges() == old num_edges() - old n.degree().
   *
   * Complexity: No more than O(num_nodes() + d*num_edges()), linear scan in valid nodes list 
   *                                                          d = n.degree() for each remove_edge operation that must be performed
   */
  size_type remove_node(const Node& n){
    if(!has_node(n)){
      return 0;
    }
    // iterate through incident edges
    std::vector<edge_proxy> n_edges = (adj_list.at(n.uid_));

    // std::cout << "Node's init degree via n_edges: " << n_edges.size() << std::endl;
    // std::cout << "Node's init degree n.degree(): " << n.degree() << std::endl;
    for(auto inci_edge = n_edges.begin(); inci_edge != n_edges.end(); ++inci_edge){
      // std::cout << "calling remove edge" << std::endl;
      remove_edge((*inci_edge).node1, (*inci_edge).node2);
      // std::cout << "Node's degree via n_edges: " << n_edges.size() << std::endl;
      // std::cout << "Node's degree n.degree(): " << n.degree() << std::endl;
    }

    // remove node from valid nodes set
    size_type node_index = n.index(); // index in valid nodes list
    valid_nodes.erase(valid_nodes_list[node_index]);

    // swap and pop in valid nodes list
    valid_nodes_list[node_index] = valid_nodes_list.back();
    valid_nodes_list.pop_back();

    // remove key from adj list
    adj_list.erase(n.uid_);

    // update num of nodes
    graph_size_ -= 1;
    assert(valid_nodes.size() == graph_size_);
    return 1;
  }

  /** remove a node and all incident edges from the graph.
   * @pre Node iterator @a n_it can be dereferenced to a valid node object n
   * @pre has_node(n) == true
   * @return a Node iterator object that points to the beginning of valid nodes
   * @post has_node(n) == false
   * @post If old has_node(n), new num_nodes() == old num_nodes() - 1.
   *       Else,               new num_nodes() == old num_nodes().
   * @post new n.degree() == 0.
   * @post new num_edges() == old num_edges() - old n.degree().
   *
   * Complexity: No more than O(num_nodes() + d*num_edges()), linear scan in valid nodes list 
   *                                                          d = n.degree() for each remove_edge operation that must be performed
   */
  node_iterator remove_node(node_iterator n_it){
    size_type index = remove_node((*n_it));
    return this->node_begin();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP