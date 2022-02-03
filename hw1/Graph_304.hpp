#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <algorithm>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
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
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->pt_map[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    /** graph_size is returned when this node is not in its original graph */
    size_type index() const {
      // HW0: YOUR CODE HERE
      for (size_type i=0; i<graph_->indices.size(); i++){
	      if (graph_->indices[i]==uid_){
	        return  i;
	      }
      }
      return graph_->indices.size()+1; //this means this node is not in the graph
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    const  node_value_type& value() const {
	    return graph_->value_map[uid_];
    }

    node_value_type& value(){
      return const_cast<node_value_type &>(static_cast<const Node &>(*this).value());
    }

    size_type degree() const{
      return (graph_->look_up[uid_]).size();
    }

    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_,0,uid_);
    }

    IncidentIterator edge_end() const{
      return IncidentIterator(graph_,degree(),uid_);
    }
    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_==n.graph_ && uid_==n.uid_;
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
      return uid_<n.uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph* graph_;
    size_type uid_;
    Node(const Graph* g,size_type uid){
      graph_=g;
      uid_=uid;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return indices.size();
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
  Node add_node(const Point& position, const node_value_type& value=node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type new_uid=pt_map.size();
    indices.push_back(new_uid);
    pt_map.push_back(position);
    validity.push_back(true);
    std::vector<size_type> edge_connect;
    look_up.push_back(edge_connect);
    value_map.push_back(value);
    return Node(this,new_uid);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    bool output=validity[n.uid_];
    //(void) n;            // Quiet compiler warning
    return output;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type node_uid=indices[i];
    //(void) i;             // Quiet compiler warning
    return Node(this,node_uid);        // Invalid node
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
      return Node(graph_,uid1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,uid2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (uid1_==e.uid1_ && uid2_==e.uid2_){
        return true;
      }
      if (uid1_==e.uid2_ && uid2_==e.uid1_){
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
      size_type less=std::min(uid1_,uid2_);
      size_type e_less=std::min(e.uid1_,e.uid2_);
      //HW0: YOUR CODE HERE
      if (less==e_less){
        size_type more=std::max(uid1_,uid2_);
        size_type e_more=std::max(e.uid1_,e.uid2_);
        return more<e_more;
      }
      return less<e_less;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    //friend class EdgeIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph* graph_;
    size_type uid1_;
    size_type uid2_;

    Edge(const Graph* g, size_type uid1, size_type uid2){
      graph_=g;
      uid1_=uid1; 
      uid2_=uid2;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_indices.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type uid1=edge_indices[i][0];
    size_type uid2=edge_indices[i][1];
    return Edge(this,uid1,uid2);        // valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::vector<size_type> a_incident=look_up[a.uid_];
    for (size_type i=0; i<a_incident.size(); i++){
      if (a_incident[i]==b.uid_){
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
    if (!has_edge(a,b)){
      std::vector<size_type> new_e;
      new_e.push_back(a.uid_);
      new_e.push_back(b.uid_);
      edge_indices.push_back(new_e);
      look_up[a.uid_].push_back(b.uid_);
      look_up[b.uid_].push_back(a.uid_);
    }
    return Edge(this,a.uid_,b.uid_);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    indices.clear();
    pt_map.clear();
    validity.clear();
    edge_indices.clear();
    look_up.clear();
    value_map.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator& operator++() {
      index++;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const{
      return index==node_iter.index && graph_==node_iter.graph_;
    }

    bool operator!=(const NodeIterator& node_iter) const{
      return !(index==node_iter.index && graph_==node_iter.graph_);
    }

    value_type operator*() const{
      return graph_->node(index);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const 

   private:
    friend class Graph;
    size_type index;
    const Graph* graph_;
    NodeIterator(const Graph* g, size_type ind): index(ind), graph_(g){}
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  NodeIterator node_begin() const{
    return NodeIterator(this, 0);
  }

  NodeIterator node_end() const{
    return NodeIterator(this, indices.size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    value_type operator*() const{
      return Edge(graph_, node_id, graph_->look_up[node_id][index]);
    }
    
    IncidentIterator& operator++(){
      index++;
      return *this;
    }

    bool operator==(const IncidentIterator& inc_iter) const{
      return node_id==inc_iter.node_id && index==inc_iter.index && graph_==inc_iter.graph_;
    }

    bool operator!=(const IncidentIterator& inc_iter) const{
      return !(node_id==inc_iter.node_id && index==inc_iter.index && graph_==inc_iter.graph_);
    }

   private:
    friend class Graph;
    friend class Node;
    size_type index;
    size_type node_id;
    const Graph* graph_;
    IncidentIterator(const Graph* g, size_type ind, size_type nid): \
                    index(ind), node_id(nid), graph_(g){}
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
    value_type operator*() const{
      return Edge(graph_,graph_->edge_indices[index][0],graph_->edge_indices[index][1]);
    }
    
    EdgeIterator& operator++(){
      index++;
      return *this;
    }

    bool operator==(const EdgeIterator& e_iter) const{
      return index==e_iter.index && graph_==e_iter.graph_;
    }

    bool operator!=(const EdgeIterator& e_iter) const{
      return !(index==e_iter.index && graph_==e_iter.graph_);
    }

   private:
    friend class Graph;
    size_type index;
    const Graph* graph_;
    EdgeIterator(const Graph* g, size_type ind): \
                    index(ind), graph_(g){}
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  EdgeIterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  EdgeIterator edge_end() const{
    return EdgeIterator(this, edge_indices.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<size_type> indices;//maps index to uid
  std::vector<Point> pt_map; //maps uid to point object
  std::vector<bool> validity;//maps uid to point validity
  std::vector<std::vector<size_type> > edge_indices;//maps index to pairs of uids
  std::vector<std::vector<size_type> > look_up;//store edge connectivity info to facilitate has_edge
  std::vector<node_value_type> value_map;//maps uid to values
};

#endif // CME212_GRAPH_HPP
