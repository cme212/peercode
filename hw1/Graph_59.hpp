#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V=int>
class Graph {
 public:
  using size_type = unsigned;
  using node_value_type = V;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  //vector to store the info and adjacency list of nodes

  struct internal_node{
    Point position_;
    size_type node_idx_;
    node_value_type val_;

    //constructor
    internal_node(Point position, size_type node_idx) : position_(position), node_idx_(node_idx){
      val_ = (node_value_type)-1;
    }
    internal_node(Point position, size_type node_idx, node_value_type val) : position_(position), node_idx_(node_idx){
      val_ = val;
    }
  };

  struct internal_edge{
    size_type node1_idx_;
    size_type node2_idx_;
    //edge value to append

    internal_edge(size_type node1_idx, size_type node2_idx): node1_idx_(node1_idx), node2_idx_(node2_idx){}

  };
  std::vector<internal_node*> nodeList_;
  std::vector<std::vector<size_type> > adjList_;

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
  Graph() : nodeList_(std::vector<internal_node*>(0)), adjList_(std::vector<std::vector<size_type> >(0)){
    // HW0: YOUR CODE HERE
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
      assert(this->graph_ != nullptr);
      assert(this->idx_ < this->graph_->nodeList_.size());
      return this->graph_->nodeList_[this->idx_]->position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx_;
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
      return this->graph_->adjList_[this->idx_].size();
    }

    // HW1: iterator part
    /** Return the start iterator for the neighbours of the node*/
    incident_iterator edge_begin() const {
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
      assert(this->graph_ == n.graph_);
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
    return this->nodeList_.size();
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    Point *newPos = new Point;
    size_type newIdx = this->size();
    Node newNode(this, newIdx);
    std::vector<size_type> newAdj;

    *newPos = position; 

    internal_node *node2Append = new internal_node(*newPos, newIdx);

    this->nodeList_.push_back(node2Append);
    this->adjList_.push_back(newAdj);
    return newNode;        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.graph_ && n.idx_ < this->size())
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
    assert(i < this->size());
    Node newNode(this, i);
    return newNode;        // Invalid node
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
      this->graph_ = graph;
      this->node1_idx_ = node1_idx;
      this->node2_idx_ = node2_idx;
    }

    Edge() {
      // empty constructor to avoid complier warning
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
      assert(this->graph_ == e.graph_);
      size_type e1n1, e1n2, e2n1, e2n2;
      e1n1 = this->node1_idx_;
      e1n2 = this->node2_idx_;

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
    const Graph* graph_;
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
    for (unsigned ii = 0;  ii < this->adjList_.size(); ii++)
    {
      for (unsigned jj = 0; jj < this->adjList_[ii].size(); jj++)
      {
        if (cnt == i)
        {
          return Edge(this, ii, this->adjList_[ii][jj]);
        }
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

    assert(node1_idx < this->size());
    assert(node2_idx < this->size());

    //check if (a, b) exsist
    for (unsigned int i = 0; i < this->adjList_[node1_idx].size(); i++)
    {
      if (adjList_[node1_idx][i] == node2_idx)
      {
        return true;
      }
    }

    //check if (b, a exist)
    for (unsigned int i = 0; i < this->adjList_[node2_idx].size(); i++)
    {
      if (adjList_[node2_idx][i] == node1_idx)
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

    assert(node1_idx < this->size());
    assert(node2_idx < this->size());

    //if not already exist, add it
    if (! this->has_edge(a, b))
    {
      this->adjList_[node1_idx].push_back(node2_idx);
      this->adjList_[node2_idx].push_back(node1_idx);
    }
    
    //just return the edge
    return Edge(this, node1_idx, node2_idx);        
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

    NodeIterator(const Graph* graph, size_type currIdx) : graph_(const_cast<Graph*>(graph)), currIdx_(currIdx){}
    

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the node denoted by this iterator*/
    Node operator*() const {
      assert(currIdx_ < graph_->size());
      Point pos = this->graph_->nodeList_[currIdx_]->position_;
      size_type nodeIdx = this->graph_->nodeList_[currIdx_]->node_idx_;
      return Node(graph_, nodeIdx);
    }

    /** Return the next node denoted by this iterator*/
    NodeIterator& operator++() {
      currIdx_ += 1;
      if (currIdx_ == graph_->size())
      {
        graph_ = nullptr;
        currIdx_ = 0;
      }
      return *this;
    }

    /** Check if the other iterator is equal to this interator */
    bool operator==(const NodeIterator& otherIter) const {
      return (graph_ == otherIter.graph_ && currIdx_ == otherIter.currIdx_);
    }

    NodeIterator& begin() {
      return NodeIterator(this->graph_, 0);
    }

    NodeIterator& end() {
      return NodeIterator(nullptr, 0);
    }

   private:
    friend class Graph;
    friend class EdgeIterator;
    // HW1 #2: YOUR CODE HERE
    Graph *graph_;
    size_type currIdx_;
  };
  


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return the start of node iterator of this graph*/
  node_iterator node_begin() const{
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
      size_type node2_idx = graph_->adjList_[node1_idx_][currIdx_];
      return Edge(graph_, node1_idx_, node2_idx);
    }

    /** Return the next edge onnected to another beighbour denoted by this iterator*/
    IncidentIterator& operator++() {
      if (currIdx_ < graph_->adjList_[node1_idx_].size()-1){
        currIdx_ += 1;
      }
      else{
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
          ++currNode_;
          if (currNode_ == graph_->node_end()){
            graph_ = nullptr;
            //currNode_ = NodeIterator(nullptr, 0);
            //currNeis_ = IncidentIterator(nullptr, 0, 0);
            return *this;
          }
          currNeis_ = IncidentIterator(graph_, currNode_.currIdx_, 0);
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
    return edge_iterator(this, 0, 0);
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
