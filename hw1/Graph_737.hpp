#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int>
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

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators which iterate over all graph edges. */
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

  }

  /** Default destructor */
  ~Graph() = default;


   // node prop struct to store additonal node information

   struct Node_Prop{
       
        private:
        friend class Graph;
        size_type degree = 0;
        std::vector<size_type> edges;
        node_value_type value = 0;

   };

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

    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->point_vec.at(id_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    node_value_type& value(){
     return graph_->node_prop[id_].value;
    }
    
    // const node_value_type& value() const;
    const node_value_type& value() const{
     return graph_->node_prop[id_].value;
    }

    // size_type degree() const;
    size_type degree() const{
     return graph_-> node_prop[id_].degree;
    }

    // incident_iterator edge_begin() const;
    incident_iterator edge_begin() const{
     incident_iterator inc_edge_iter;
     inc_edge_iter.id_=id_;
     inc_edge_iter.idx=0;
     inc_edge_iter.graph_=graph_;
     return inc_edge_iter;
    }

    // incident_iterator edge_end() const;
    incident_iterator edge_end() const{
     incident_iterator inc_edge_iter;
     inc_edge_iter.id_=id_;
     inc_edge_iter.idx=graph_->node_prop[id_].edges.size();
     inc_edge_iter.graph_=graph_;
     return inc_edge_iter;
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (id_==n.id_){
       if(&(graph_->point_vec[id_])==&(n.graph_->point_vec[id_])){
        return true;
       }
       return false; 
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
      if(id_<n.id_){
       return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer to graph class that contains this node
    Graph* graph_;
    // Node id
    size_type id_;
    // bool describing wether this node has been validated
    bool valid = 0;


    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // store the point
    point_vec.emplace_back(position);
    // initialize invalid point
    Node node_;
    // alter node attributes
    node_.graph_=this;
    node_.id_=num_nodes();
    node_.valid=true;
    // add node to graph
    node_vec.emplace_back(node_);
    node_prop.emplace_back(Node_Prop());
    return node_; 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check id
    if (n.id_<num_nodes()){
     if (((n.graph_->point_vec[n.id_].x) == (point_vec[n.id_].x)) and
         ((n.graph_->point_vec[n.id_].y) == (point_vec[n.id_].y)) and
         ((n.graph_->point_vec[n.id_].z) == (point_vec[n.id_].z))){
      return true;
     }
     return false;
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
    assert(i<num_nodes());
    return node_vec[i];
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      size_type id_node1=node1id;
      return graph_->node_vec[id_node1];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type id_node2=node2id;
      return graph_->node_vec[id_node2];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
       size_type id_node1=node1id;
       size_type id_node2=node2id;

       if (graph_->node_vec[id_node1]==graph_->node_vec[e.node1id]){
        if (graph_->node_vec[id_node2]==graph_->node_vec[e.node1id])
         return true;
       }else if (graph_->node_vec[id_node1]==graph_->node_vec[e.node2id]){      
        if (graph_->node_vec[id_node2]==graph_->node_vec[e.node1id])
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
      if (id_<e.id_){
       return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // edge nodes
    size_type node1id;
    size_type node2id;

    // edge id
    size_type id_;

    // give access to graph
    Graph* graph_;

    // Validity
    bool valid = 0;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edge_vec[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // create set with node ids
    std::set<size_type> node_set={a.id_,b.id_};
    if(edge_map.count(node_set))
     return true;
    else
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
    std::set<size_type> node_set={a.id_,b.id_};
    if (has_edge(a,b)){
     return edge_vec[ edge_map[node_set]];
    }
    else{
     // initialize invalid edge
     Edge edge_;
     // set attributes
     edge_.id_=num_edges();
     edge_.graph_=this;
     edge_.node1id=a.id_;
     edge_.node2id=b.id_; 
     node_prop[a.id_].degree += 1;
     node_prop[b.id_].degree += 1;
     node_prop[a.id_].edges.emplace_back(edge_.id_);
     node_prop[b.id_].edges.emplace_back(edge_.id_);
     edge_.valid=1;
     edge_vec.emplace_back(edge_);
     edge_map[node_set]=edge_.id_;
     return edge_;   
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_vec.clear();
    edge_vec.clear();
    edge_map.clear();
    point_vec.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    Node operator*() const {
     return graph_->node_vec[idx];
    }
    // NodeIterator& operator++()
    NodeIterator& operator++(){
     idx++;
     return *this;
    }
    // bool operator==(const NodeIterator&) const
    bool operator==(const NodeIterator& node_iter) const{
     return idx==node_iter.idx;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type idx;
    const Graph* graph_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  node_iterator node_begin() const{
   node_iterator node_;
   node_.idx=0;
   node_.graph_= this;
   return node_;
  }
  // node_iterator node_end() const
  node_iterator node_end() const{
   node_iterator node_;
   node_.idx=node_vec.size();
   node_.graph_= this;
   return node_;
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
    // Edge operator*() const
    Edge operator*() const{
     return graph_->edge_vec[graph_->node_prop[id_].edges[idx]];
    }
    // IncidentIterator& operator++()
    IncidentIterator& operator++(){
     idx++;
     return *this;
    }
    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& inc_edge_iter) const{
     return (id_==inc_edge_iter.id_ && idx==inc_edge_iter.idx);
    }
   
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type id_;
    size_type idx;
    const Graph* graph_;
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
    // Edge operator*() const
    Edge operator*() const{
     return graph_->edge_vec[idx];
    }
    // EdgeIterator& operator++()
    EdgeIterator& operator++(){
     idx++;
     return *this;
    }
    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& edge_iter) const{
     return idx==edge_iter.idx;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type idx;
    const Graph* graph_; 
 };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const{
   edge_iterator edge_;
   edge_.idx=0;
   edge_.graph_= this;
   return edge_;
  }
  // edge_iterator edge_end() const
  edge_iterator edge_end() const{
   edge_iterator edge_;
   edge_.idx=edge_vec.size();
   edge_.graph_= this;
   return edge_;
  }

 private:
  // Store nodes
  std::vector<node_type> node_vec;

  // Store node information

  std::vector<Node_Prop> node_prop;
 
  // Store edges
  std::vector<edge_type> edge_vec;

  // Store edge nodes
  std::map<std::set<size_type>,size_type> edge_map;

  // Store points
  std::vector<Point> point_vec;

};

#endif // CME212_GRAPH_HPP
