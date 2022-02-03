#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct Node_attr;
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Type of node value. */
  using node_value_type=V;
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
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
    clear(); // free up all dynamically allocated memory
  }

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
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(G_!=nullptr);// invalid node does not have a position
      // retrieve position from node_vec which stores node attributes
      return G_->node_vec[node_idx]->node_pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // Calling on invalid node may return undefined node_idx
      return node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Retrieve the value associated with the node from its attribute
    node_value_type& value(){
      assert(G_!=nullptr);
      return G_->node_vec[node_idx]->node_val;
    }

    // Retrieve the value associated with the node from its attribute
    // Return as const to avoid modification of the value
    const node_value_type& value() const{
      assert(G_!=nullptr);
      return G_->node_vec[node_idx]->node_val;
    }

    // Return the size of incident edge vectors stored in node attribute
    size_type degree() const {
      if(G_==nullptr)
        return 0; // Invalid Node have zero degree
      return G_->node_vec[node_idx]->incid_edge.size();
    }
    
    // Start the incident iterator from node indexed as 0
    incident_iterator edge_begin() const {
      return IncidentIterator(G_, G_->node_vec[node_idx], 0);
    }
    
    // Incident iterator end should point one past over all nodes
    incident_iterator edge_end() const {
      return IncidentIterator(G_, G_->node_vec[node_idx], 
G_->node_vec[node_idx]->incid_edge.size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      /** Check whether they share the same Graph pointer
         to the graph where they belongs to and whether
         the two nodes have the same index. */
      if((G_==nullptr)&&(n.G_==nullptr))
        return true; // two invalid nodes will be considered equal
      if((G_==n.G_) && (node_idx==n.node_idx)) {
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
      // order based on node index
      // make sure to compare within same graph or between invalid nodes
      if (*this==n)
        return false;
      if (G_==n.G_ && node_idx<n.node_idx) {
        return true;
      }
      return false;
    }


    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type node_idx=size_type(); // node index
    Graph* G_=nullptr; // Pointer to the graph where the node belongs to

    /** Private node constructor to construct a valid Node */
    Node(const Graph* G, size_type idx) {
      G_=const_cast<Graph*>(G);
      node_idx=idx;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_vec.size();
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
  Node add_node(const Point& position, const node_value_type& val= node_value_type()) {
    // HW0: YOUR CODE HERE
    Node_attr* n=new Node_attr(position,val,size());
    node_vec.push_back(n);// insert new node attribute pointer
    return Node(this,size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // index of node attr pointer in the node_vec is equal to the node_idx
    // of the node
    // the node should exist in this graph
    if (n.G_==this&&n.node_idx<size())
      return true;
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
    // if query an existing node, return a copy of Node with such input index
    if (i<size()) {
      return Node(this,i);
    }
    return Node(); // return invalid node if not exist
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // Check whether two edges are in the same graph and
      // share the same connected nodes
      if (G_==nullptr && e.G_==nullptr)
        return true; // consider both invalid edges as equivalent
      if (G_!=e.G_) { // valid edges should have a non null graph pointer
        return false;
      }else if ((*n1==*e.n1 && *n2==*e.n2)||(*n1==*e.n2 && *n2==*e.n1)) {
        return true;
      }else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // make sure to compare within the same graph, or between invalid edges
      if (*this==e)
        return false;
      if (G_==e.G_ && edge_idx<e.edge_idx) { // order by edge index
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // pointers to two nodes connected by the edge
    const Node* n1;
    const Node* n2;
    size_type edge_idx=size_type(); // index of edge
    Graph* G_=nullptr; // Graph pointer to the graph where the edge belongs to

    /** Private constructor to construct a valid edge. */
    Edge(Graph* G, const Node& a, const Node& b, size_type idx) {
      G_=const_cast<Graph*>(G);
      n1=&a;
      n2=&b;
      edge_idx=idx;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i<num_edges()) {
      return *(edge_vec[i]); // return deferenced edge pointer
    }
    return Edge(); // return an invalid Edge if cannot find input indexed edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    /* use helper function which returns a pointer to such edge cononecting
       the two nodes, or nullptr if no such edge. */
    Edge* temp=get_edge(a,b);
    if (temp!=nullptr)
      return true;
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
    // make sure input valid nodes
    Edge* temp=get_edge(a,b);
    if (temp!=nullptr) { // return pointer to current edge if already exists
      return Edge(this,a,b,temp->edge_idx);
    }

    size_type edge_count=num_edges();

    Edge* new_edge=new Edge(this, a, b, edge_count);
    // Store edge index that is incident to the two nodes
    node_vec[a.node_idx]->incid_edge.push_back(edge_count);
    node_vec[b.node_idx]->incid_edge.push_back(edge_count);
    edge_vec.push_back(new_edge);
    node_2_edge.insert({std::make_pair(a.node_idx, b.node_idx),
    new_edge->edge_idx});
    return *new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // free up dynamically allocated memory first
    for (auto& v: node_vec) {
      if (v!=nullptr) {
        delete v;
      }
    }
    for (auto& v: edge_vec) {
      if (v!=nullptr) {
        delete v;
      }
    }
    // clear up containers
    node_vec.clear();
    node_2_edge.clear();
    edge_vec.clear();

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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Deference a NodeIterator. */
    Node operator*() const{
      assert(G_!=nullptr); // should not deference an invalid node iterator
      return G_->node(iter_idx);
    }

    /** Forward increment the NodeIterator by incrementing iterator index. */
    NodeIterator& operator++(){
      ++iter_idx;
      return *this;
    }

    /** Equivalent node iterators should belong to the same graph and
        point to the node with the same node index. */
    bool operator==(const NodeIterator& iter) const{
      if(G_==iter.G_ && iter_idx==iter.iter_idx){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type iter_idx=size_type();
    Graph* G_=nullptr;
    // Valid private node iterator constructor
    NodeIterator(const Graph* G, size_type idx){
      G_=const_cast<Graph*>(G);
      iter_idx=idx;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // Node iterator should begin to point to node indexed 0
  node_iterator node_begin() const{
    return  NodeIterator(this,0);
  }

  /** Update iterator index to be one past the total number
      of nodes in the graph. */
  node_iterator node_end() const{
    return NodeIterator(this,size());
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
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** The edge obtained by deferencing the incident iterator should have the
        the node that spawns the iterator as node1. */
    Edge operator*() const {
      /* Retrieve the existing edge incident to the node where the iterator
         points to. */
      assert(G_!=nullptr&&n_attr!=nullptr); // should not deref invalid incident iterator
      Edge temp=G_->edge((*n_attr).incid_edge[outedge_idx]);
      /* If node1 of the existing edge is not the node that spawns the
         iterator, return an edge with flipped node1 and node2. */ 
      if(temp.node1()==G_->node((*n_attr).node_idx)){
        return temp;
      }else{
        Edge temp2=G_->add_edge(temp.node2(),temp.node1());
        return temp2;
      }
    }

    /* Increment the incident iterator by moving onto the next edge
       incident to the node pointed to by the iterator. */
    IncidentIterator& operator++(){
      outedge_idx++;
      return *this;
    }
    
    /* Equivalent incident iterators should belong to the same graph,
       point to the same node, and currently examine the same edge incident
       to that node. */
    bool operator==(const IncidentIterator& iter) const {
      if(G_==iter.G_ && n_attr==iter.n_attr &&
outedge_idx==iter.outedge_idx){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Node_attr* n_attr=nullptr;
    Graph* G_=nullptr;
    /** index of edge w.r.t out_edge vector from node attribute that
      the IncidentIterator is currently point to. */
    size_type outedge_idx=size_type();

    // Private valid incident iterator constructor
    IncidentIterator(const Graph* G,const Node_attr* attr, size_type idx){
      G_=const_cast<Graph*>(G);
      n_attr=const_cast<Node_attr*>(attr);
      outedge_idx=idx;
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
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Dereference operator that returns the edge pointed to by the iterator
    Edge operator*() const {
      assert(G_!=nullptr); // should not dereference an invalid edge iterator
      return G_->edge(edge_idx);
    }

    /** Forward increment operation on EdgeIterator, move onto the next edge
        in the graph with incremented edge index. */
    EdgeIterator& operator++() {
      edge_idx++; 
      return *this;
    }

    /** Equivalent edge iterators should belong to the same graph and point
        to the edge with the same edge index. */
    bool operator==(const EdgeIterator& iter) const {
      if(G_==iter.G_ && edge_idx==iter.edge_idx){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* G_=nullptr;
    size_type edge_idx=size_type();

    // Valid private edge iterator constructor
    EdgeIterator(const Graph* G, size_type idx) {
      G_=const_cast<Graph*>(G);
      edge_idx=idx;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // Edge iterator begins to point to the first edge in the graph
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  // Edge iterator ends one past over all edges in the graph
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct Node_attr {
      Point node_pos;
      node_value_type node_val;
      size_type node_idx;
      // store indices of edges incident to the node
      std::vector<size_type> incid_edge;

      Node_attr(const Point& P, const node_value_type& val, size_type idx) {
        node_pos=Point(P.x, P.y, P.z);
        node_val=val;
        node_idx=idx;
      }
  };
  std::vector<Node_attr*> node_vec; // vector of Node pointers
  // map with the pair of node indices as key and store edge index
  std::map<std::pair <size_type,size_type>,size_type> node_2_edge;
  // vector of edge pointer indexed as edge idx
  std::vector<Edge*> edge_vec;

  /** Get a pointer to edge with node1==a and node2==b or vice
     versa in the graph.
   * @pre @a a and @a b are valid nodes of the graph
     @return the pointer to edge if for some @a i, edge(@a i)
      connects @a a and @a b.*/
  Edge* get_edge(const Node& a, const Node& b) const {
    // make sure input valid nodes in the same graph
    assert(a.G_==this&&b.G_==this);
    // search whether node indices key pair exist
    auto search1=node_2_edge.find(std::make_pair(a.node_idx,b.node_idx));
    auto search2=node_2_edge.find(std::make_pair(b.node_idx,a.node_idx));
     if (search1!=node_2_edge.end()) {
      return edge_vec[search1->second];
    } else if (search2!=node_2_edge.end()) {
      return edge_vec[search2->second];
    }
    return nullptr ; // return nullptr if not found
  }
};

#endif // CME212_GRAPH_HPP

