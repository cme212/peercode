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


template <typename V>
class Graph {

    
    



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
    
  using node_value_type = V;
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
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): size_(0), next_uid_(0){
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
    Node()
      {
    }
      
     const node_value_type& value () const {return value_;}
     void distance(int di) {this->dist = di;}
     
      
     int getDistance() {
         return dist;
     };
      
      
      
     size_type degree() const {
         return adjEdges.size();
     };
      
      void add_adjEdge(Edge& e){
          adjEdges.push_back(e);
      }
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
        
      if (uid_ > g_->size()){
          assert(false);}
      return g_->nodes[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    incident_iterator edge_begin() const {
        return incident_iterator(const_cast<std::vector<Edge>*>(&adjEdges), 0);
    };
      
      
    incident_iterator edge_end() const {
        return incident_iterator(const_cast<std::vector<Edge>*>(&adjEdges), adjEdges.size());
    };
        
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      return (n.g_ == g_) && (n.uid_ == uid_);
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
      return n.uid_ < uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    Graph* g_;  
    size_type uid_;
    const node_value_type& value_;
//     size_type degree_;
    std::vector<Edge> adjEdges;
    int dist;
      
    
      
    Node(const Graph* g, size_type uid, const node_value_type& val = node_value_type()): g_(const_cast<Graph*>(g)), uid_(uid), value_(val), dist(0){}
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
      
      
    // HW0: YOUR CODE HERE
    return next_uid_;
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
  const Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    nodes.push_back(position);
    size_ ++;
    next_uid_ ++;
    vertices.push_back(Node(this, next_uid_-1, val));
    return Node(this, next_uid_-1, val);
      
      
      // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
 
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return this == n.g_;
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
//     return Node(this, i);        // Invalid node
      return (vertices[i]);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return node1_;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return node2_;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
        
        
      //HW0: YOUR CODE HERE
      return ((e.node1_ == node1_) && (e.node2_ == node2_)) || ((e.node1_ == node2_) && (e.node2_ == node1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return e.eid_ < eid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
      Node& node1_;
      Node& node2_;
      unsigned int eid_;
      Graph* graph_;  
    friend class Graph;
      Edge(const Graph*graph, size_type eid, Node& n1, Node& n2):graph_(const_cast<Graph*>(graph)), eid_(eid), node1_(n1), node2_(n2){}
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };
//     Edge(const Graph*graph, size_type eid, Node n1, Node n2):graph_(const_cast<Graph*>(graph)), eid_(eid), node1_(n1), node2_(n2){}
   
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
//     (void) i;             // Quiet compiler warning
    return edges[i];        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
//     (void) a; (void) b;   // Quiet compiler warning
    for (Edge e: edges){
        bool d1 = (e.node1_ == a) && (e.node2_ == b);
        bool d2 = (e.node1_ == b) && (e.node2_ == a);
            
        if (d1 || d2)
            return true;
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
  void set_dist(int ind, int d1){
      vertices[ind].distance(d1);
  }
    
    
 void print_dist(){
     for (auto ni = node_begin(); ni != node_end(); ++ni){
      
         Node v = *ni;
         std::cout<<v.getDistance()<<"distance"<<v.index()<<"index"<<std::endl;
     }
 }
  Edge add_edge(Node& a, Node& b) {
      
    for (Edge e: edges){
        bool d1 = (e.node1_ == a) && (e.node2_ == b);
        bool d2 = (e.node1_ == b) && (e.node2_ == a);
            
        if (d1 || d2)
            return e;
        
    }
      
      Edge e = Edge(this, edges.size(), a, b);
      Edge e2 = Edge(this, edges.size(), b,a);

      int id1 = a.index();
      int id2 = b.index();
      a.add_adjEdge(e);
      b.add_adjEdge(e2);
      vertices[id1].add_adjEdge(e);
      vertices[id2].add_adjEdge(e2);
              
      edges.push_back(e);
    // HW0: YOUR CODE HERE
//     (void) a, (void) b;   // Quiet compiler warning
    return e;        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      nodes.clear();
      edges.clear();
    // HW0: YOUR CODE HERE
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
      
    
    NodeIterator(): nodes_vector(nullptr), curIndex(0){
        
    }
      
    NodeIterator(std::vector<Node>* vec, int idx): nodes_vector(vec), curIndex(idx){}
      
    Node operator*() const {
        return (*nodes_vector)[curIndex];
    }
      
    NodeIterator& operator++(){
        curIndex ++;
        return *this;
    }
      
    
    bool operator==(const NodeIterator& node_iteror) const {
        
        
        return (node_iteror.curIndex == curIndex);
    }
      
    bool operator!=(const NodeIterator& node_iteror) const {
        return !((*this) == node_iteror);
    
    }
      
    int get_size(){
        return (*nodes_vector).size();
        
    
    }
      
    int get_index(){
        return curIndex;}
      
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Node;
    friend class Graph;
    std::vector<Node>* nodes_vector;
    int curIndex;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
      
      return node_iterator(const_cast<std::vector<Node>*> (&vertices), 0);
  }
  node_iterator node_end() const  {
      return node_iterator(const_cast<std::vector<Node>*> (&vertices), vertices.size());
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
    IncidentIterator():incidence_e(nullptr), edge_ind(0){
    }
      
    IncidentIterator(std::vector<Edge>* eds, int idx):incidence_e(eds), edge_ind(idx){}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
        return (*incidence_e)[edge_ind];
        
    }
    
    IncidentIterator& operator++(){
        
        edge_ind ++;
        return *this;
    }
      
    bool operator==(const IncidentIterator& iteror) const{
        return (edge_ind == iteror.edge_ind) ;
    
    }
      
    bool operator !=(const IncidentIterator& iteror) const {
        return !((*this) == iteror);
    }
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    friend class Node;
    friend class Edge;
    std::vector<Edge>* incidence_e;
    int edge_ind;
  
    
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
    EdgeIterator(): n1_(NodeIterator()), i1_(IncidentIterator()), g_(nullptr){
        
    }
    EdgeIterator(NodeIterator n1,  const Graph* g) {
        this->n1_ = n1;
        this->g_ = g;
        this->i1_ = IncidentIterator();
        if (this->n1_ != g->node_end()){ 

            this->i1_ = (*n1_).edge_begin();
            
        }

    }
    EdgeIterator(NodeIterator& n1, IncidentIterator& i1, const Graph* g): n1_(n1), i1_(i1), g_(const_cast<Graph*>(g)){}

    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
        return *i1_;
    }
      
    EdgeIterator& operator++(){
       
            bool invalid = true;
            while (invalid) {
                if (i1_ != (*n1_).edge_end()){
                    std::cout<<"increment edge"<<std::endl;
                    ++i1_;
                    if(i1_ == (*n1_).edge_end()) continue;
                }
                
                else{  
                    std::cout<<"increment node";
                    ++n1_;
                    if (n1_ == g_->node_end()){
                        std::cout<<"end of the graph nodes"<<std::endl;
                        invalid = false;
                        break;
                    }
                    else{
                        std::cout<<(*n1_).index()<<"node_index"<<std::endl;
                        i1_ = (*n1_).edge_begin();
                        std::cout<<"go to another node"<<(*n1_).index()<<std::endl;
                        
                    }
                }
                
                if (!invalid) {break;}
                else{
                    if ((*n1_).edge_begin() != (*n1_).edge_end()){
                        
//                         std::cout<<(*n1_).index()<<"n1"<< (*i1_).node2().index()<<"n2"<<std::endl;
//                         std::cout<<((*n1_).index() < (*i1_).node2().index())<<"check valid or not"<<std::endl;
                        
                        if ((*n1_).index() < (*i1_).node2().index()){
                            std::cout<<"find a valid edge"<<std::endl;
                            invalid = false;
                            break;
                        }
                    }   
                } 
            }
        
        std::cout<<"iter"<<std::endl;
        
        return (*this);
    }
    bool operator==(const EdgeIterator& eiter) const {
        if (n1_ == g_->node_end() && eiter.n1_ == g_->node_end()) return true;
        return (n1_ == eiter.n1_) && (i1_ == eiter.i1_) && (g_ == eiter.g_);
    }
      
    bool operator !=(const EdgeIterator& eiter) const {
        
        return !((*this) == eiter);
    }

   private:
    friend class Graph;
    NodeIterator n1_;
    IncidentIterator i1_;
//     Graph* g_;
    const Graph* g_;
    
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
      NodeIterator start_Node = node_begin();
      return EdgeIterator(start_Node, this);
  
  
  }
  edge_iterator edge_end() const {
  
      NodeIterator end_Node = node_end();
      
//       return EdgeIterator();
      return EdgeIterator(end_Node, this);
  
  }

 private:
    
    
    std::vector<Point> nodes;
    std::vector<Edge> edges;
    std::vector<Node> vertices;
    unsigned int size_;
    unsigned int next_uid_;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
