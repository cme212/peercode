#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @ class Graph
 * @ brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_value_type = E;

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
  /** struct to store the elements to a node.
   * idx: index of node
   * point: position of node
   * value: auxiliary value of node.
   */
  struct NodeElement{
   size_type idx;
   Point point;
   node_value_type value;
  };


  /** struct to store the elements to a node.
   * idx: index of edge
   * node1_idx: uid of node1 of current edge
   * node2_idx: uid of node2 of current edge
   * value: auxiliary value of edge. 
   */
  struct EdgeElement{
   size_type idx;
   size_type node1_idx;
   size_type node2_idx;
   edge_value_type value;
  };

  /** Below is the representation R(G = <N, E>). 
   * nodes and edges are vectors to store the information.
   * node_adj maps each node to a set of edges incidenting to it.
   * The use of set facilitate the usage of search, insert 
   * and delete operations. 
   * idx2uid is seperating uid and idx of node.
   */
  std::vector<NodeElement> nodes;
  std::vector<EdgeElement> edges;
  std::unordered_map<size_type, std::set<size_type>> node_adj;
  std::vector<size_type> idx2uid;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    /** Do nothing here to construct a empty graph. */
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
  class Node: private totally_ordered<Node> {
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
      // Doesn't do anything for the invalid constructor.
    }
    
    
    /** Non const function of positions. */
    Point& position(){
      assert(idx_ < graph_->size());
      size_type uid = graph_ -> idx2uid[idx_];
      return graph_ -> nodes[uid].point;
    }



    /** Return this node's position. */
    const Point& position() const {
      // The node needs to be valid.
      assert(idx_ < graph_->size());
      size_type uid = graph_ -> idx2uid[idx_];
      return graph_ -> nodes[uid].point;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }


    // HW1: YOUR CODE HERE
    /** Return this node's value.*/
    node_value_type &value(){
      assert(idx_ < graph_->size());
      size_type uid = graph_ -> idx2uid[idx_];
      return graph_ -> nodes[uid].value;
    }

    /** Return this node's value. */
    const node_value_type& value() const{
      assert(idx_ < graph_->size());
      size_type uid = graph_ -> idx2uid[idx_];
      return graph_ -> nodes[uid].value;
    }


    /**Set the value
     */
    void setValue(node_value_type new_value){
      assert(idx_ < graph_->size());
      size_type uid = graph_ -> idx2uid[idx_];
      graph_ -> nodes[uid].value = new_value;
    }
    

    // Supply definitions AND SPECIFICATIONS for:
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    /**@brief Degree of node is defined as the number of
     * edges that connects to the current node.
     *@return The number of edges that incident to the current node.
     */
    size_type degree() const{
      size_type uid = graph_ -> idx2uid[idx_];
      return graph_ -> node_adj[uid].size();
    }


    /**@brief Get the begin iterator for the incident iterator.
     * In my data strucuture, each node maps a list for the edges.
     * The edge list is traversed and the first element is the starting
     * element and the last one is the end element.
     *@return An incident iterator that refers to the begining
     * value. 
     */
    incident_iterator edge_begin() const{
      size_type uid = graph_ -> idx2uid[idx_];
      return IncidentIterator(graph_, 
                              idx_, 
                              graph_->node_adj[uid].begin());
    }


    /**@brief get the end iterator for the incident iterator.
     *@return An incident iterator that refers to the end value.
     */
    incident_iterator edge_end() const{
      size_type uid = graph_ -> idx2uid[idx_];
      return IncidentIterator(graph_, 
                              idx_, 
                              graph_->node_adj[uid].end());
    }
    


    /** Test whether this node and _n_ are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->idx_ == n.idx_ && this->graph_ == n.graph_;
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
      return (this->graph_ < n.graph_) or 
             (this->graph_ == n.graph_ and this->idx_ < n.idx_);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph *graph_;
    size_type idx_;

    // Constructor only available in graph class.
    Node(const Graph* graph_input, const size_type id_input)
         : graph_(const_cast<Graph*>(graph_input)), idx_(id_input){}
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return idx2uid.size();
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
                const node_value_type & v = node_value_type()) {
    size_type num = size();
    NodeElement node_added = {num, position, v};
    nodes.push_back(node_added);
    Node new_node = Node(this, num);
    idx2uid.push_back(num);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ and n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, i);
  }

  
  /** Remove the current node if it exists. If node exists, all 
   * incident edges need to be removed as well.
   * @param[in] _node_ the node to be removed from the graph.
   * @return 1 if the node exist and remove the node from the graph. 
   *         The returned value is 0 if the node doesn't exist.
   * If node doesn't exists, return 0 and no other changes happen.
   * If node exists, @post size() = size() - 1. Edges with current 
   * nodes needs to be removed from the map.
   * 
   * Complexity: O(log (node.degree()))
   */
  size_type remove_node(const Node& node){
    if(!has_node(node)){
      return 0;
    }
    size_type curr = node.index();
    size_type removed_node_uid = idx2uid[curr];
    idx2uid[curr] = idx2uid.back();
    nodes[idx2uid[curr]].idx = curr;
    idx2uid.pop_back();
    
    for(auto i:node_adj[removed_node_uid]){
      size_type uid_1 = edges[i].node1_idx;
      size_type uid_2 = edges[i].node2_idx;
      node_adj[uid_1].erase(i);
      node_adj[uid_2].erase(i);
    }

    return 1;
  }

  /** Remove the current node if it exists. If node exists, all 
   * incident edges need to be removed as well.
   * @param[in] _n_it_ An node iterator.
   * @return A valid node iterator as the node.begin().
   * If node doesn't exists, no other changes happen and return begin iterator.
   * If node exists, @post size() = size() - 1. Edges with current 
   * nodes needs to be removed from the map.
   * 
   * Complexity: O(log (node.degree()))
   */
  node_iterator remove_node(node_iterator n_it){
    remove_node((*n_it));
    return node_begin();
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      /* Doesn't do anything for an invalid edge. */
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, (graph_->nodes[node1_idx_]).idx); 
    }

    double length() const{
      double length = norm(graph_->nodes[node1_idx_].point -
                      graph_ -> nodes[node2_idx_].point);
      return length;
    }

    /** Return this edge's value.*/
    edge_value_type &value(){
      return graph_ -> edges[idx_].value;
    }

    /** Return this edge's value. */
    const edge_value_type& value() const{
      return graph_ -> edges[idx_].value;
    }


    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, (graph_->nodes[node2_idx_]).idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->graph_ == e.graph_ and idx_ == e.idx_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->graph_ < e.graph_) or (this->graph_ == e.graph_
              and idx_ < e.idx_);
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_idx_;
    size_type node2_idx_;

    // We define an id for edge like the id for node.
    size_type idx_;

    // Constructor only available for graph class. 
    Edge(const Graph* graph_input, const size_type node1_id_input, 
      const size_type node2_id_input, const size_type id_input):
        graph_(const_cast<Graph*>(graph_input)),
        node1_idx_(node1_id_input), 
        node2_idx_(node2_id_input), 
        idx_(id_input){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges(){
    size_type count = 0;
    for(size_type i : idx2uid){
      count += node_adj[i].size();
    }
    return size_type(count/2);
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < edges.size() and i>=0);
    return Edge(this, edges[i].node1_idx,
                edges[i].node2_idx,
                edges[i].idx);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O((log num_edges())), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    // find_edge is another function defined below.
    // It is used to find the id of edge that connects two nodes.
    // The complexity is O(num_edges()).
    size_type i = find_edge(a, b);
    if(i < edges.size())
      return true;
    return false;
  }



  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre _a_ and _b_ are distinct valid nodes of this graph
   * @return an Edge object e with _e_._node1_() == _a_ and _e_._node2_() == _b_
   * @post has_edge(_a_, _b_) == true
   * @post If old has_edge(_a_, _b_), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(_i_) might not
   * equal to new edge(_i_). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type & e = edge_value_type()) {
    assert(has_node(a) && has_node(b) && ! (a==b));
    size_type num = edges.size();
    size_type uid_a = idx2uid[a.index()];
    size_type uid_b = idx2uid[b.index()];
    size_type i = find_edge(a, b);
    if(i == num){
      Edge new_edge = Edge(this, uid_a, uid_b, num);
      EdgeElement edge_added = {num, uid_a, uid_b, e};
      edges.push_back(edge_added);
      node_adj[uid_a].insert(i);
      node_adj[uid_b].insert(i);
      return new_edge;
    }else{
      return Edge(this, uid_a, uid_b, num);
    }
  }



  /** Given two node and find the id of the edge that connects them.
   * This is a helpful function and
   * @pre _node_a_ and  _node_b_ are different valid nodes
   * @return the id of edge if it exists. Otherwise, return the number of edges.
   * 
   * Complexity: O(log (num_edges())).
   */
  size_type find_edge(const Node& node_a, const Node& node_b) const{
    // We need to ensure that node_a/b are different and valid.
    assert(has_node(node_a) && has_node(node_b));
    if(node_a == node_b){
      return edges.size();
    }
    size_type uid_a = idx2uid[node_a.index()];
    size_type uid_b = idx2uid[node_b.index()];
    if(node_adj.find(uid_a) != node_adj.end()){
      for(size_type i:node_adj.at(uid_a)){
        if(edges[i].node1_idx == uid_b ||
           edges[i].node2_idx == uid_b)
          return i;
      }
    }
    return edges.size();
  }



  /** Remove the current edge if it exists. 
   * @param[in] _e_ the edge to be removed from the graph.
   * @return 1 if the edge exist and remove the edge from the graph. 
   *         The returned value is 0 if the edge doesn't exist.
   * If edge doesn't exists, return 0 and no other changes happen.
   * If edge exists, @post num_edges() = num_edges() - 1. 
   * 
   * Complexity: O(log (num_edges()))
   */
  size_type remove_edge(const Edge &e){
    Node n1 = node(nodes[e.node1_idx].idx);
    Node n2 = node(nodes[e.node2_idx].idx);
    if(!has_edge(n1, n2)){
      return 0;
    }
    size_type curr_node1_uid = e.node1_idx;
    size_type curr_node2_uid = e.node2_idx;
    node_adj[curr_node1_uid].erase(e.idx_);
    node_adj[curr_node2_uid].erase(e.idx_);
    return 1;
  }



  /** Remove the current edge if it exists. 
   * @param[in] _n1_, _n2_ the nodes of edge
   * to be removed from the graph.
   * @return 1 if the edge exist and remove the edge from the graph. 
   *         The returned value is 0 if the edge doesn't exist.
   * If edge doesn't exists, return 0 and no other changes happen.
   * If edge exists, @post num_edges() = num_edges() - 1. 
   * 
   * Complexity: O(log (num_edges()))
   */
  size_type remove_edge(const Node& n1, const Node & n2){
    if(!has_edge(n1, n2)){
      return 0;
    }
    size_type curr = find_edge(n1, n2);
    size_type uid_1 = idx2uid[n1.index()];
    size_type uid_2 = idx2uid[n2.index()];
    node_adj[uid_1].erase(curr);
    node_adj[uid_2].erase(curr);
    return 1;
  }


  /** Remove the current edge if it exists. 
   * @param[in] _e_it_ the iterator of edge
   * to be removed from the graph.
   * @return edge_begin().
   * If edge exists, @post num_edges() = num_edges() - 1. 
   * 
   * Complexity: O(log (num_edges()))
   */
  edge_iterator remove_edge(const edge_iterator e_it){
    remove_edge((*e_it));
    return edge_begin();
  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_adj.clear();
    idx2uid.clear();
  }
  
  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator>{

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

    /** @brief Get the value that the current iterator refers to.
     * @return A node object. 
     * Complexity: O(1)
     */
    Node operator*() const{
      return Node(graph_, idx_);
    }
    

    /**@brief Increment the current iterator.
     *@return A iterator with different value to reference.
     *complexity: O(1)
     */
    NodeIterator &operator++(){
      idx_ = idx_ + 1;
      return (*this);
    }
 

    /**@brief Check if two iterator are the same.
     *@param[in] The iterator (*this) and it. 
     *@return return true if the two are the same, otherwise, 
     * return false.
     *Complexity: O(1)
     */
    bool operator==(const NodeIterator & it) const{
      return this->graph_ == it.graph_ && this->idx_ == it.idx_;
    }

   private:
    friend class Graph;
    Graph *graph_;
    size_type idx_;

    // Constructor for the graph and id.
    NodeIterator(const Graph *graph_input, const size_type id_input):
                 graph_(const_cast<Graph*>(graph_input)), idx_(id_input){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /**@brief Define the start iterator.
   *@return return a node that serves as the start iterator.
   *complexity:O(1)
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  


  /**@brief Define the end iterator.
   *@return return a node that serves as the end iterator.
   *complexity:O(1)
   */
  node_iterator node_end() const{
    return NodeIterator(this, this->size());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
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

    /** @brief Get the value that the current iterator refers to.
     * @return An edge object. 
     * @post The node1 of the returned edge must be equal to the 
     * id_node_ below.
     * Complexity: O(1)
     */
    Edge operator*() const{
      size_type id_curr_edge = (*edge_iter_);
      EdgeElement e = graph_ -> edges[id_curr_edge];
      if(e.node1_idx == id_node_) 
      return Edge(graph_, e.node1_idx, e.node2_idx,  id_curr_edge);
      else return Edge(graph_, e.node2_idx, e.node1_idx, id_curr_edge);
    }

    /**@brief Increment the iterator.
     *@return An incident iterator object with incremented iterator.
     * complexity: O(1)
     */
    IncidentIterator &operator++(){
      edge_iter_++;
      return (*this);
    }


    /**@brief Check if the two incident iterators are the same.
     *@return True if they are the same or false.
     */
    bool operator==(const IncidentIterator& it) const{
      return this->graph_ == it.graph_ && this->id_node_ == it.id_node_
             && this->edge_iter_ == it.edge_iter_;
    }


   private:
    friend class Graph;
    Graph *graph_;
    // Current node that edges incident to.
    size_type id_node_;
    std::set<size_type>::iterator edge_iter_;
    IncidentIterator(const Graph* graph, const size_type id_node, 
                    std::set<size_type>::iterator edge_iter):
                    graph_(const_cast<Graph*>(graph)), 
                    id_node_(id_node), 
                    edge_iter_(edge_iter){}
  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
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
    /**@brief Get the current edge
     *@return Get the edge that the current edge iterator refers to
     * Complexity: O(1)
     */
    Edge operator*() const{
      return *(i_it_);
    }


    /**@brief Increment the current edge iterator.
     *@return Get the edge iterator that refers to the next value.
     * Complexity: O(1)
     */
    EdgeIterator &operator ++(){
      increment();
      return (*this);
    }
    

    /**@brief Check if two Edge iterators are the same.
     *@return True if the same or false.
     */
    bool operator == (const EdgeIterator & it) const{
      return this -> graph_ == it.graph_ 
             && this->n_it_ == it.n_it_
             && this->i_it_ == it.i_it_;
    }


   private:
    friend class Graph;
    Graph *graph_;
    node_iterator n_it_;
    incident_iterator i_it_;
    // Constructor
    EdgeIterator(const Graph* graph, 
                 const node_iterator& n_it, 
                 const incident_iterator &i_it)
                : graph_(const_cast<Graph*>(graph)), 
                  n_it_(n_it), i_it_(i_it){}
    void increment(){
      while(1){
        ++i_it_;
        if(i_it_ == (*n_it_).edge_end()){
          ++n_it_;
          if(n_it_ == graph_->node_end()){
            i_it_ = (*(graph_->node_begin())).edge_begin();
            break;
          }
          i_it_ = (*n_it_).edge_begin();
        }
        if((*i_it_).node1_idx_ + (*i_it_).node2_idx_ 
            < 2*graph_->idx2uid[(*n_it_).idx_]){
          break;
        }
      }
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /**@brief In the data structure, edges are stored with unique id
   * in a list. The first element in the list is used as the starting
   * element and the last one is the end element.
   * @return Starting edge iterator.
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const{
    node_iterator curr_n_it = node_begin();
    incident_iterator curr_i_it = (*curr_n_it).edge_begin();
    return EdgeIterator(this, curr_n_it, curr_i_it);
  }


  /**brief Return the end edge iterator.
   *@return End edge iterator.
   *Complexity: O(1)
   */
  edge_iterator edge_end() const{
    node_iterator curr_n_it = node_end();
    incident_iterator curr_i_it = (*node_begin()).edge_begin();
    return EdgeIterator(this, curr_n_it, curr_i_it);
  }
};

#endif // CME212_GRAPH_HPP
