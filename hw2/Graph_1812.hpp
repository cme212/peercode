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
template <typename V=int, typename E=double>
class Graph {
 public:
  using size_type = unsigned;
  using node_value_type = V;

  using edge_len_type = double;
  using edge_value_type = E;

 private:
  //friend class Edge;

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  //vector to store the info and adjacency list of nodes

  struct internal_node{
    Point position_;
    size_type id4user_; // node_idx_ should be id4user_. idx is it's index inthe vector
    node_value_type val_;

    //constructor
    internal_node(Point position, size_type id4user) : position_(position), id4user_(id4user){
      val_ = node_value_type {};
    }
    internal_node(Point position, size_type id4user, node_value_type val) : position_(position), id4user_(id4user){
      val_ = val;
    }
  };

  struct internal_edge{
    size_type node1_idx_;
    size_type node2_idx_;
    //edge value to append

    internal_edge(size_type node1_idx, size_type node2_idx): node1_idx_(node1_idx), node2_idx_(node2_idx){}

  };

  std::vector<internal_node*> nodeList_; //index by idx which is unique and not shown to user
  std::vector<std::vector<std::pair<size_type, edge_value_type> > > adjList_; //this should be indexed by idx!

  //edge iterator & node iterator must be re-written! along with add node & add edge
  std::vector<size_type> idu2idx_; //index by id_for_user, return the corresponding idx



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
  

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodeList_(std::vector<internal_node*>(0)), adjList_(std::vector<std::vector<std::pair<size_type, edge_value_type> > >(0)){
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
    clear(); //release the resource gained by "new"
  };

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
    Node() {} //construct invalid node

    Node(const Graph* graph, unsigned idx): graph_(graph), idx_(idx) {
      this->graph_ = graph;
      this->idx_ = idx;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_ != nullptr);
      assert(idx_ < graph_->nodeList_.size());
      return graph_->nodeList_[idx_]->position_;
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      assert(graph_ != nullptr);
      assert(idx_ < graph_->nodeList_.size());
      return graph_->nodeList_[idx_]->position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return this->idx_;
      assert(graph_->nodeList_[idx_]->id4user_ < graph_->size()); //or it will be invalid node.
      return graph_->nodeList_[idx_]->id4user_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the value of this node, which can be modified*/
    node_value_type& value() {
      return this->graph_->nodeList_[this->idx_]->val_;
    }

    /** Return the value of this node, which cannot be modified*/
    const node_value_type& value() const{
      return this->graph_->nodeList_[this->idx_]->val_;
    }

    /** Return the degree of this node*/
    size_type degree() const{
      return this->graph_->adjList_[idx_].size();
    }

    // HW1: iterator part
    /** Return the start iterator for the neighbours of the node*/
    incident_iterator edge_begin() const {
      if (degree() == 0){
        return IncidentIterator(nullptr, 0, 0);
      }
      return IncidentIterator(graph_, idx_, 0);
    }

    /** Return the end iterator for the neighbours of the node*/
    incident_iterator edge_end() const{
      return IncidentIterator(nullptr, 0, 0);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->graph_ == n.graph_ && this->idx_ == n.idx_)
      {
        return true;
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
      // HW0: YOUR CODE HERE
      if (this->graph_ != n.graph_){
        return graph_ < n.graph_;
      }
      return this->idx_ < n.idx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph* graph_;
    size_type idx_;
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //return this->nodeList_.size();
    return idu2idx_.size();
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
  Node add_node(const Point& position, node_value_type val = node_value_type()) {
    // HW0: YOUR CODE HERE
    //Point *newPos = new Point;
    size_type newIdu = size();
    size_type newIdx = nodeList_.size();

    Node newNode(this, newIdx); //node should always initialized by idx
    std::vector<std::pair<size_type, edge_value_type>> newAdj;

    //*newPos = position; 

    internal_node *node2Append = new internal_node(position, newIdu, val);//internal_node should record the idu
    nodeList_.push_back(node2Append);
    adjList_.push_back(newAdj);
    
    idu2idx_.push_back(newIdx); //record the currnode's idx. the id4user is the index of idu2idx.
    return newNode;        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.graph_ && n.idx_ < nodeList_.size() && nodeList_[n.idx_]->id4user_ < size())
    {
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
    // HW0: YOUR CODE HERE
    assert(i < size());
    Node newNode(this, idu2idx_[i]);
    return newNode;        // Invalid node
  }

  /**
   * @brief Given a node @n n, remove it from the graph. The node is invalidated after removal.
   * 
   * @param[in] @n n a Node to be removed 
   * @return True if _n_ was in the graph and is successfully removed, else false
   * @post has_node(@n n) == false
   * @post If old has_node(@n n), new num_nodes() == old num_edges() - 1.
   *       Else,                  new num_nodes() == old num_edges().
   * 
   * Complexity: O(num_edges + num_nodes)
   */

  size_type remove_node (const Node& n){
    //deleting all the incident edges 
    if (!has_node(n)){
      return false;
    }
    //incident_iterator nei = n.edge_begin();
    //while (nei != n.edge_end()){
    //  nei = remove_edge(nei);

    while (n.degree() > 0){
      Edge nei = *(n.edge_begin());
      remove_edge(nei);
    }
  
    size_type nodeIdx = n.idx_;
    size_type nodeIdu = nodeList_[nodeIdx]->id4user_;

    size_type lastNodeIdu = size()-1;
    size_type lastNodeIdx = idu2idx_[lastNodeIdu];

    //deactivating the node to be removed
    nodeList_[nodeIdx]->id4user_ = -1; //which will be a very large number if we use unsigned

    //let the last node move to the place of the node to be removed
    nodeList_[lastNodeIdx]->id4user_ = nodeIdu;
    idu2idx_[nodeIdu] = lastNodeIdx;

    idu2idx_.pop_back();

    return true;
  }

  /**
   * @brief Given a node iterator n_it, remove the node it points to from the graph
   * 
   * @param[in] @n n a Node to be removed 
   * @return True if _n_ was in the graph and is successfully removed, else false
   * @post has_node(@n n) == false
   * @post If old has_node(@n *n_it), new num_nodes() == old num_edges() - 1.
   *       Else,                      new num_nodes() == old num_edges().
   * 
   * Complexity: O(num_edges + num_nodes)
   */

  node_iterator remove_node (node_iterator n_it){
    assert(n_it != node_end());

    size_type flag = remove_node(*n_it);
    if (flag == false || size() == 0 || n_it.currIdu_ >= size()){
      return NodeIterator(nullptr, 0);
    }

    return n_it;
    
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
    Edge(const Graph* graph, size_type node1_idx, size_type node2_idx) {
      // HW0: YOUR CODE HERE
      graph_ = const_cast<Graph*>(graph);
      node1_idx_ = node1_idx;
      node2_idx_ = node2_idx;
    }

    Edge() {
      // empty constructor to avoid complier warning
    }

    
    edge_value_type& value(){
      for (size_type j=0; j< graph_->adjList_[node1_idx_].size(); j++){
        if (graph_->adjList_[node1_idx_][j].first == node2_idx_){
          return (*graph_).adjList_[node1_idx_][j].second;
        }
      }

      //should never reach here
      exit(1);
    }

    void setValue(edge_value_type val2set){
      assert(graph_ != NULL);
      for (size_type j=0; j< graph_->adjList_[node1_idx_].size(); j++){
        if (graph_->adjList_[node1_idx_][j].first == node2_idx_){
          (*graph_).adjList_[node1_idx_][j].second = val2set;
          //graph_->adjList_[node1_idx_][j].second = val2set;
          //std::cout << graph_->adjList_[node1_idx_][j].first <<std::endl;
        }
      }
    }

    edge_value_type value_noref(){
      for (size_type j=0; j< graph_->adjList_[node1_idx_].size(); j++){
        if (graph_->adjList_[node1_idx_][j].first == node2_idx_){
          return (*graph_).adjList_[node1_idx_][j].second;
        }
      }

      return 0;
    }

    const edge_value_type& value() const{
      for (size_type j=0; j < graph_->adjList_[node1_idx_].size(); j++){
        if (graph_->adjList_[node1_idx_][j].first == node2_idx_){
          return graph_->adjList_[node1_idx_][j].second;
        }
      }

      return 0;
    }
    

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node1_idx_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node2_idx_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (this->graph_ == e.graph_ && ((this->node1_idx_ == e.node1_idx_ && this->node2_idx_ == e.node2_idx_) \
         || (this->node1_idx_ == e.node2_idx_ && this->node2_idx_ == e.node1_idx_)))
         {
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
      //HW0: YOUR CODE HERE

      //first compare the smaller node_id, if equivalent, campare the greater node_id
      //assert(this->graph_ == e.graph_);
      if (graph_ != e.graph_){
        return graph_ < e.graph_;
      }

      size_type e1n1, e1n2, e2n1, e2n2;
      e1n1 = node1_idx_;
      e1n2 = node2_idx_;

      if (e1n2 < e1n1)
      {
        size_type tmp = e1n1;
        e1n1 = e1n2;
        e1n2 = tmp;
      }

      e2n1 = e.node1_idx_;
      e2n2 = e.node2_idx_;

      if (e2n2 < e2n1)
      {
        size_type tmp = e2n1;
        e2n1 = e2n2;
        e2n2 = tmp;
      }

      if (e1n1 < e2n1)
      {
        return true;
      }

      if (e1n1 == e2n1 && e1n2 < e2n2)
      {
        return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class EdgeIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_idx_;
    size_type node2_idx_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type cnt = (size_type)0;
    for (unsigned int i = 0; i < this->adjList_.size(); i++)
    {
      cnt += this->adjList_[i].size();
      
    }
    return cnt/2;
  }

  
  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < this->num_edges());
    size_type cnt = (size_type)0;
    for (unsigned ii = 0;  ii < adjList_.size(); ii++)
    {
      for (unsigned jj = 0; jj < adjList_[ii].size(); jj++)
      {
        if (cnt == i)
        {
          return Edge(this, ii, adjList_[ii][jj].first);
        }
        cnt++;
      }
    }
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == b.graph_);

    size_type node1_idx = a.idx_;
    size_type node2_idx = b.idx_;

    assert(has_node(a));
    assert(has_node(b));

    //check if (a, b) exsist
    for (unsigned int i = 0; i < this->adjList_[node1_idx].size(); i++)
    {
      if (adjList_[node1_idx][i].first == node2_idx)
      {
        return true;
      }
    }

    //check if (b, a exist)
    for (unsigned int i = 0; i < this->adjList_[node2_idx].size(); i++)
    {
      if (adjList_[node2_idx][i].first == node1_idx)
      {
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == b.graph_);

    size_type node1_idx = a.idx_;
    size_type node2_idx = b.idx_;

    assert(has_node(a));
    assert(has_node(b));

    //if not already exist, add it
    if (! this->has_edge(a, b))
    {
      std::pair<size_type, edge_value_type> pair1(node2_idx, edge_value_type());
      std::pair<size_type, edge_value_type> pair2(node1_idx, edge_value_type());
      adjList_[node1_idx].push_back(pair1);
      adjList_[node2_idx].push_back(pair2);
      //adjList_[node2_idx][adjList_[node2_idx].size()-1].second = 1;
      //adjList_[node1_idx][adjList_[node1_idx].size()-1].second = 1;
    }
    
    //just return the edge
    return Edge(this, node1_idx, node2_idx);        
  }

  /**
   * @brief Given two nodes, remove the edge between the nodes from the graph if it exists
   *        After deleting, the edge is invalidated.
   * @param[in,out] node1 a node in one end of the edge to be removed
   * @param[in,out] node2 another node in the other end of the edge to be removed
   * @return True if the edge was in the graph and is deleted successfully
   * 
   * @pre node1 and node2 are valid nodes in the graph
   * @post if has_edge(@node1 node1, @node2 node2), new num_edges() == old num_edges() - 1.
   *       Else,                                    new num_edges() == old num_edges()
   * @post has_edge(@node1 node1, @node2 node2) is false.
   * Complexity: O(d), d is the maximun degree of the node in the graph.
   */

  size_type remove_edge(const Node& node1, const Node& node2){
    if (!has_node(node1) || !has_node(node2) || !has_edge(node1, node2)){
      return false;
    }

    size_type node1Idx = node1.idx_;
    size_type node2Idx = node2.idx_;

    size_type len = adjList_[node1Idx].size();
    for (size_type i = 0; i < len; i++)
    {
      if (adjList_[node1Idx][i].first == node2Idx)
      {
        adjList_[node1Idx][i].first = adjList_[node1Idx][len-1].first;
        adjList_[node1Idx][i].second = adjList_[node1Idx][len-1].second;
        adjList_[node1Idx].pop_back();
        break;
      }
    }

    len = adjList_[node2Idx].size();
    for (size_type i = 0; i < len; i++)
    {
      if (adjList_[node2Idx][i].first == node1Idx)
      {
        adjList_[node2Idx][i].first = adjList_[node2Idx][len-1].first;
        adjList_[node2Idx][i].second = adjList_[node2Idx][len-1].second;
        adjList_[node2Idx].pop_back();
        break;
      }
    }
    return true; //delete witout error
  }

  /**
   * @brief Given an edge @e e, remove the edge from the graph if it exists.
   *        After deleting, the edge is invalidated.
   * @param[in] @e _e_ a node in one end of the edge to be removed
   * @return True if the edge was in the graph and is deleted successfully
   * 
   * 
   * @post if has_edge(@node1 e.node1(), @node2 e.node2()), new num_edges() == old num_edges() - 1.
   *       Else,                                            new num_edges() == old num_edges()
   * @post has_edge(@node1 e.node1(), @node2 e.node2()) is false.
   * Complexity: O(d), d is the maximun degree of the node in the graph.
   */

  size_type remove_edge (const Edge& e){
    node_type node1 = e.node1();
    node_type node2 = e.node2();

    return remove_edge(node1, node2);
  }
  
  /**
   * @brief Given an edge interator, remove the edge it points to from the graph if it exists.
   *        After deleting, the edge iterator is invalidated.
   * @param[in,out] e_it a node in one end of the edge to be removed
   * @return The invalidated edge iterator
   * 
   * @post if has_edge(@node1 (*e_it).node1(), @node2 (*e_it).node2()), new num_edges() == old num_edges() - 1.
   *       Else,                                            new num_edges() == old num_edges()
   * @post has_edge(@node1 (*e_it).node1(), @node2 (*e_it).node2()) is false.
   * Complexity: O(d), d is the maximun degree of the node in the graph.
   */
  
  edge_iterator remove_edge (edge_iterator e_it){
    assert(e_it != edge_end());
    remove_edge(*e_it);
    if (adjList_[e_it.node1_idx_].size() == 0)
    {
      return EdgeIterator(nullptr, 0, 0);
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
    this->nodeList_.clear();
    for (unsigned i = 0; i < this->adjList_.size(); i++)
    {
      this->adjList_[i].clear();
    }
    this->adjList_.clear();

    this->idu2idx_.clear();
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


    NodeIterator() {
    }

    NodeIterator(const Graph* graph, size_type currIdu) : graph_(const_cast<Graph*>(graph)), currIdu_(currIdu){}
    

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the node denoted by this iterator*/
    Node operator*() const {
      
      //assert(currIdx_ < graph_->size());
      //Point pos = this->graph_->nodeList_[currIdx_]->position_;
      //size_type nodeIdx = this->graph_->nodeList_[currIdx_]->node_idx_; //nodeIdx is actually the currIdx, when not considering the removal of node
      //return Node(graph_, nodeIdx);

      assert(currIdu_ < graph_->size());
      size_type nodeIdx = graph_->idu2idx_[currIdu_];
      return Node(graph_, nodeIdx);
    }

    /** Return the next node denoted by this iterator*/
    NodeIterator& operator++() {
      //assert(currIdx_ < graph_->size());
      currIdu_ += 1;
      if (graph_ == nullptr || currIdu_ >= graph_->size())
      {
        graph_ = nullptr;
        currIdu_ = 0;
      }
      return *this;
    }

    /** Check if the other iterator is equal to this interator */
    bool operator==(const NodeIterator& otherIter) const {
      return (graph_ == otherIter.graph_ && currIdu_ == otherIter.currIdu_);
    }

    /*
    NodeIterator& begin() {
      return NodeIterator(this->graph_, 0);
    }

    NodeIterator& end() {
      return NodeIterator(nullptr, 0);
    }
    */

   private:
    friend class Graph;
    friend class EdgeIterator;
    // HW1 #2: YOUR CODE HERE
    Graph *graph_;
    size_type currIdu_;
  };
  


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return the start of node iterator of this graph*/
  node_iterator node_begin() const{
    if (size() == 0){
      return NodeIterator(nullptr, 0);
    }
    return NodeIterator(this, 0);
  }

  /** Return the end of node iterator of this graph*/
  node_iterator node_end() const{
    return NodeIterator(nullptr, 0);
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
    IncidentIterator(const Graph* graph, size_type node1_idx, size_type currIdx = 0){
      graph_ = const_cast<Graph*>(graph);
      node1_idx_ = node1_idx;
      currIdx_ = currIdx;
    } 

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge connected to a beighbour denoted by this iterator*/
    Edge operator*() const {
      assert(currIdx_ < graph_->adjList_[node1_idx_].size());
      size_type node2_idx = graph_->adjList_[node1_idx_][currIdx_].first;
      return Edge(graph_, node1_idx_, node2_idx);
    }

    /** Return the next edge onnected to another beighbour denoted by this iterator*/
    IncidentIterator& operator++() {
      currIdx_ += 1;

      if (graph_ == nullptr || currIdx_ >= graph_->adjList_[node1_idx_].size()){
        graph_ = nullptr;
        node1_idx_ = 0;
        currIdx_ = 0;
      }
      return *this;
    }

    /** Check if the two itertor are equal*/
    bool operator==(const IncidentIterator& otherIter) const {
      return (graph_ == otherIter.graph_ && node1_idx_ == otherIter.node1_idx_ && currIdx_ == otherIter.currIdx_);
    }

   private:
   // HW1 #3: YOUR CODE HERE
    friend class Graph;
    friend class EdgeIterator;
    Graph *graph_;
    size_type node1_idx_;
    size_type currIdx_;
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

    EdgeIterator(const Graph* graph, size_type currNodeIdx, size_type currNeisIdx) {
      graph_ = const_cast<Graph*>(graph);
      currNode_ = NodeIterator(graph, currNodeIdx);
      currNeis_ = IncidentIterator(graph, currNodeIdx, currNeisIdx);
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge in the graph denoted by this iterator*/
    Edge operator*() const {
      return *currNeis_;
    }

    /** Return the next edge in the graph denoted by this iterator*/
    EdgeIterator& operator++() {
      do{
        ++currNeis_;
        if (currNeis_.graph_ == nullptr){
          do{
            ++currNode_;
          }while(currNode_ != graph_->node_end() && (*currNode_).degree() == 0);
          
          if (currNode_ == graph_->node_end()){
            graph_ = nullptr;
            //currNode_ = NodeIterator(nullptr, 0);
            //currNeis_ = IncidentIterator(nullptr, 0, 0);
            return *this;
          }
          currNeis_ = IncidentIterator(graph_, graph_->idu2idx_[currNode_.currIdu_], 0);
        }     
      } while ((*currNeis_).node1_idx_ > (*currNeis_).node2_idx_);
      return *this;
    }

    /** Check if the two iterator are equal*/
    bool operator==(const EdgeIterator& otherIter) const {
      return (graph_ == otherIter.graph_ && currNode_ == otherIter.currNode_ && currNeis_ == otherIter.currNeis_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph *graph_;
    NodeIterator currNode_;
    IncidentIterator currNeis_;
    
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Return the start of edge interator of this graph*/
  edge_iterator edge_begin() const {
    for (size_type i = 0; i < size(); i++)
    {
      if (node(i).degree() != 0)
      {
        return edge_iterator(this, i, 0);
      }
    }
    return edge_iterator(nullptr, 0, 0);
  }

  /** Return the end of edge interator of this graph*/
  edge_iterator edge_end() const {
    return edge_iterator(nullptr, 0, 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
