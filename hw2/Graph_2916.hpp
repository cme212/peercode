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

template <typename V, typename E>
class Graph : private totally_ordered<Graph<V,E>> {
 typedef V node_value_type; 
 typedef E edge_value_type;

 private:

 public:
  std::vector<Point> nodePoints;
  std::vector<V> val; //For use with NodeData type

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

  std::map<std::vector<size_type>, size_type> Edges;
  std::vector<size_type> vec;  //vec[i] = uid of node at index i
  //(namely, uid of node corresponding to nodePoints[i])
  std::vector<size_type> uidInd; //uidInd[i] = index of node with
  //uid = i

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  friend class NodeIterator;
  /** Construct an empty graph. */
  Graph(): size_(0) {
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
    Graph* graph_;
    size_type index_;
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

    //@pre caller is a valid node in the graph
    //@return the value associated with the node
    node_value_type& value() {
      return (*graph_).val[index_];
    }

    //@pre caller is a valid node in the graph
    //@return the value associated with the node
    const node_value_type& value() const {
      return graph_->val[index_];
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodePoints[index_];
    }

    /** Return this node's position. */
    Point& position() {
      return graph_->nodePoints[index_];
    }

    //@return the degree of the node (number of incident edges)
    size_type degree() const {
      size_type cumDeg = (*graph_).IncEdges[index_].size();
      return cumDeg;
    }

    //@return iterator pointing at first incident edge
    incident_iterator edge_begin() const {   
      return incident_iterator(graph_, (*graph_).IncEdges[index_], 
        index_, 0);
    }

    //@return iterator pointing at last incident edge
    incident_iterator edge_end() const {
      return incident_iterator(graph_, (*graph_).IncEdges[index_], 
        index_, (*this).degree());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return (*this).graph_->uidInd[index_];
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
      return(index_ == n.index_ && graph_ == n.graph_);
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
      return(index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class nodeColor;

    Node(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {} 
    //index_ is the uid of the node
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
  Node add_node(const Point& position, const V& v) {
    nodePoints.push_back(position);
    size_ += 1;
    vec.push_back(size_-1);
    uidInd.push_back(size_-1);
    Node addedNode = Node(this, size_ - 1); 
    ((*this).val).push_back(v);
    return addedNode;
  }

 /** Removes a node from the graph if it exists
   * @param[in] n the node to be removed
   * @return 1 if the node exists in current graph, 0 otherwise
   * 
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - degree(n)
   * @post node n passed in is invalidated if it exists in graph
   *
   * Complexity: O(num_nodes()) operations.
   */
  size_type remove_node(const Node& n) {
    //Check to see if node exists, returns 0 immediately if not
    size_type ind = 0;
    for (size_type i=0; i<size_; i++) {
      if (vec[i] == n.index_){
        ind = 1;
      }
    }
    if (ind == 0) {return 0;}

    //First remove incident edges
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto e = *ei;
      remove_edge(e);
    }
    
    //Now remove the node
    for (size_type i=0; i<size_; i++) {
      if (vec[i] == n.index_) {
        //Update internal mapping of uid's to indices and 
        //indices to uid's
        auto j = vec[size_-1];
        vec[i]=j;
        uidInd[j]=i;

        //Update points and values so that the ith entry 
        //corresponds to the new node 
        nodePoints[i] = nodePoints[size_-1];
        val[i] = val[size_-1];
      }
    }
    nodePoints.pop_back();
    //val.pop_back();
    size_ -= 1;
    return 1;
  }

  /** Removes a node from the graph if it exists.
   * @param[in] n_it  valid iterator to node to remove
   * @return iterator to first node if size_>0, otherwise
   * the end iterator
   * 
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - degree(n)
   * @post n_it is invalid 
   *
   * Complexity: O(num_nodes()) operations.
   */
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    auto ind = remove_node(n);
    if (size_ > 0) {
      return this->node_begin();
    } else {
      return this->node_end();
    }
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index_ <= size_ - 1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    auto j = (*this).vec[i];
    //std::cout << j << std::endl;
    return Node(this, j);
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

    size_type index1_;
    size_type index2_;

    Edge(const Graph* graph, size_type index1, size_type index2) 
    : graph_(const_cast<Graph*>(graph)), index1_(index1), index2_(index2) {}
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index2_);     // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool a = (index1_ == e.index1_) && (index2_ == e.index2_);
      bool b = (index1_ == e.index2_) && (index2_ == e.index1_);
      //bool c = 
      return ((a || b) && (graph_ == e.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //First compare smallest indices, and only if they are equal
      //compare the largest

      //std::cout << index1_ << index2_ << std::endl; 
      //std::cout << e.index1_ << e.index2_ << std::endl; 
      if (graph_ != e.graph_) {
        return true;
      }
      else if (index1_ > e.index1_) {
        return false;  
      } else if (index1_ < e.index1_) {
        return true;
      } else  {
        return (index2_ < e.index2_);
      }
    }

    //@pre caller is a valid edge in the graph
    //@return the value associated with the edge
    edge_value_type& value() {
      std::vector<size_type> search{node1().index_, node2().index_};
      return (graph_->edgeVals.find(search))->second;
    }

    //@pre caller is a valid edge in the graph
    //@return the value associated with the edge
    const edge_value_type& value() const {
      std::vector<size_type> search{node1().index_, node2().index_};
      return (graph_->edgeVals.find(search))->second;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;    
    Graph* graph_;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return Edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
   */
  Edge edge(size_type i) const {
    edge_type desiredEdge;
    for (auto& it : Edges) {
        if (it.second == i+1) {
          desiredEdge = Edge(this, it.first[0], it.first[1]);
          break;
        }
    }
    return desiredEdge;
  }

 /** Return a 2 entry vector whose 1st entry is the minimum.
   * @pre 2 unsigned integers a and b
   *
   */
  std::vector<size_type> minimum(size_type a, size_type b) const {
    std::vector<size_type> ordIndices;
    if (a < b) {
      ordIndices.push_back(a);
      ordIndices.push_back(b);
    }
    else { 
      ordIndices.push_back(b);
      ordIndices.push_back(a);
    }
    return ordIndices;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(log(num_edges()))
   */
  bool has_edge(const Node& a, const Node& b) const {
    //To ensure we only search for the edge (a,b) once, edges
    //are stored with nodes ordered with increasing index
    std::vector<size_type> ordered = minimum(a.index_, b.index_);
    auto min = ordered[0];
    auto other = ordered[1];
    std::vector<size_type> search{min, other};
    return (Edges.find(search) != Edges.end());
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
   * Complexity: O(log(num_edges)))
   */
  Edge add_edge(const Node& a, const Node& b, const E& e) {
    //E e;
    if (IncEdges.find(a.index_) == IncEdges.end()) {
      IncEdges[a.index_].push_back(b.index_);
    } else {
      std::vector<size_type> adj = IncEdges[a.index_];
      bool ind =false;
      for (size_type i=0; i< adj.size(); i++) {
        if (adj[i] == b.index_) {
          ind = true;
        }
      }
      if (ind == false) {
        IncEdges[a.index_].push_back(b.index_);
      }
    }

    if (IncEdges.find(b.index_) == IncEdges.end()) {
      IncEdges[b.index_].push_back(a.index_);
    } else {
      std::vector<size_type> adj = IncEdges[b.index_];
      bool ind =false;
      for (size_type i=0; i< adj.size(); i++) {
        if (adj[i] == a.index_) {
          ind = true;
        }
      }
      if (ind == false) {
        IncEdges[b.index_].push_back(a.index_);
      }
    }

    if (has_edge(a,b)) {
      return(Edge(this, a.index_, b.index_));
    } else {
      std::vector<size_type> ordered = minimum(a.index_, b.index_);
      size_type min = ordered[0];
      size_type other = ordered[1];
      std::vector<size_type> search1{min, other};
      std::vector<size_type> search2{other, min};
      //In our map data structure, keys are node pairs corresponding to edges,
      //values are number of edges added
      Edges[search1] = this->num_edges() + 1;
      edgeVals[search1] = e;
      edgeVals[search2] = e;
      return Edge(this, a.index_, b.index_);
    }

  }

 /** Removes an edge from the graph if it exists
   * @param[in] n1, n2 endpoint nodes of edge to be removed
   * @return 1 if the edge exists in current graph, 0 otherwise
   * 
   * @post If old has_edge(n1, n2), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * @post has_edge(e.n1(), e.n2()) == false
   *
   * Complexity: O(num_nodes + num_edges()) operations.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1, n2)){
      return 0;
    }

    std::vector<size_type> adj1 = IncEdges[n1.index_];
    std::vector<size_type> adj2 = IncEdges[n2.index_]; 
    //Remove the "other node" from the Indicent Edges data structure
    for(auto curr=adj1.begin(); curr != adj1.end(); ++curr) {
      if (*curr == n2.index_) {
        adj1.erase(curr);  
        break;
      }
    }
    for(auto curr=adj2.begin(); curr != adj2.end(); ++curr) {
      if (*curr == n1.index_) {
        adj2.erase(curr);
        break;
      }
    }

    //Now removes the edge from the main Edges data structure
    std::vector<size_type> ordered = minimum(n1.index_, n2.index_);
    for(auto it = Edges.begin(); it != Edges.end(); ++it) {
        if(it->first == ordered) {
            Edges.erase(it); 
            break;
        } 
    }

    //Finally removes any edgevalues associated with the edges
    size_type first = ordered[0];
    size_type second = ordered[1];
    std::vector<size_type> other{second, first};
    for(auto it = edgeVals.begin(); it != edgeVals.end(); ++it) {
        if(it->first == ordered) {
            edgeVals.erase(it);  //'it' is invalidated so need to 
            //start a new loop for the other edge orientation
            break;
        } 
    }
    for(auto it = edgeVals.begin(); it != edgeVals.end(); ++it) {
        if(it->first == other) {
            edgeVals.erase(it); 
            break;
        } 
    }

    return 1;
  }

  /** Removes an edge from the graph if it exists
   * @param[in] e edge to be removed
   * @return 1 if the edge exists in current graph, 0 otherwise
   * 
   * @post If old has_edge(n1, n2), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * @post has_edge(e.node1(), e.node2()) == false
   *
   * Complexity: O(num_nodes + num_edges()) operations.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Removes an edge from the graph if it exists
   * @param[in] e_it valid iterator to edge to be removed
   * @return the beginning edge_iterator if there are still edges
   * after removal, otherwise the ending one
   * 
   * @pre valid edge_iterator 
   * @post If old has_edge(n1, n2), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * @post e_it passed in is invalidated if its edge exists in graph
   *
   * Complexity: O(num_nodes + num_edges()) operations.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    size_type a = remove_edge(e);
    if (a && Edges.size() > 0) {
      return edge_iterator(this, 0);
    } else if (Edges.size() == 0) {
      return (*this).edge_end();
    } else {
      return e_it;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodePoints.clear();
    size_ = 0;
    Edges.clear();
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
    const Graph* graph_;
    int index_;
    NodeIterator(const Graph* graph, int index) : 
    graph_(graph), index_(index) {}

    //@return the current node pointed to by iterator
    Node operator*() const {
        return Node(graph_, index_);
    }

    //@return a reference to an iterator pointed at the next node
    //@post index_ is one larger
    NodeIterator& operator++() {
        index_++;
        return *this;
    }

    //@param nodeIt  The iterator being compared to
    //@return true if both iterators pointing at same node, false 
    //otherwise
    bool operator==(const NodeIterator& nodeIt) const {
        return(index_ == nodeIt.index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /*Starts the node iterating process.
  @return An iterator to the first (index 0) node*/ 
  node_iterator node_begin() const {
      return node_iterator(this, 0);
  }

  /*Marks the end of the node iterating process.
  @return An iterator to the one-after-last (index nodePoints.size()-1) node*/ 
  node_iterator node_end() const {
      return node_iterator(this, nodePoints.size());
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

    const Graph* graph_;
    std::vector<size_type> otherInd_;
    size_type baseIndex_;
    size_type index_;
    IncidentIterator(const Graph* graph, std::vector<size_type> otherInd,
     int baseIndex, int index) : 
    graph_(graph), otherInd_(otherInd), 
    baseIndex_(baseIndex), index_(index) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    //@return The current edge
    Edge operator*() const {
        return Edge(graph_, baseIndex_, otherInd_[index_]); 
      }
      
    /*@return iterator pointng at next incident edge
    @post index_ is one larger (namely at next adjacent node)*/
    IncidentIterator& operator++() {
      index_++;
      return *this;
    }

    /*@param incIt  The incident iterator being compared to
    @return true if both edges are the same, false otherwise*/
    bool operator==(const IncidentIterator& incIt) const {
      return(baseIndex_ == incIt.baseIndex_ &&
        //index_ == otherInd_.size());
        index_ == incIt.index_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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


    const Graph* graph_;
    size_type index_;
    EdgeIterator(const Graph* graph, size_type index) : 
    graph_(graph), index_(index) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    //@return The current edge
    Edge operator*() const {
      return((*graph_).edge(index_));
    }

    /*@return iterator pointng at next edge in sequence
    @post index_ is one larger (namely at next edge)*/
    EdgeIterator& operator++() {
      index_++;
      return *this;
    }

    /*@param edgeIt  The edgeiterator being compared to
    @return true if both edges are the same, false otherwise*/
    bool operator==(const EdgeIterator& edgeIt) const {
      return (index_ == edgeIt.index_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  //@return iterator pointing at first edge 
  edge_iterator edge_begin() {
    return EdgeIterator(this, 0);
  }

  //@return iterator pointing at last edge   
  edge_iterator edge_end() {
    return EdgeIterator(this, Edges.size());
  }

 private:

  size_type size_;
  std::map<std::vector<size_type>, E> edgeVals;
  std::map<size_type, std::vector<size_type>> IncEdges;

};

#endif // CME212_GRAPH_HPP
