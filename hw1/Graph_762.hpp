#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <map>
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
class Graph{
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;

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

  private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  //attribute: vector of nodes
  // std::vector<std::tuple<Point*,node_value_type>> nodes_;
  
  //Using maps to reduce runtime
  std::map<size_type,std::tuple<Point*,node_value_type>> nodes_;
  //attribute: edge list of vectors of node indices of each edge
  
  // std::vector<std::vector<size_type>> edge_list_;
  //Also using maps to reduce runtime
  //map[a] will contain all nodes that a is connected to
  std::map<size_type,std::vector<size_type>> edge_list_;

  public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  Graph():nodes_(std::map<size_type,std::tuple<Point*,node_value_type>>()),
          edge_list_(std::map<size_type,std::vector<size_type>>()){}

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
  class Node: private totally_ordered <Node> {
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
      //return dereferenced Point
      // return *(graph->nodes_[std::node_idx]);
      return *(std::get<0>(graph->nodes_[node_idx]));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      return node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      return std::get<1>(graph->nodes_.at(node_idx));
    };
    //set value function to assign value in shortest_path
    void set_value(node_value_type& new_value){
      std::get<1>(graph->nodes_[node_idx]) = new_value;
    };
    const node_value_type& value() const{
      return (std::get<1>(graph->nodes_[node_idx]));
    };
    size_type degree() const{
      //return degree of this node
      return graph->edge_list_[node_idx].size();
    };
    incident_iterator edge_begin() const{
      return IncidentIterator(graph,node_idx,0);
    };
    incident_iterator edge_end() const{
      return IncidentIterator(graph,node_idx,graph->edge_list_.at(node_idx).size());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //If two nodes have same graph and same index then they should be equal
      if (this->graph == n.graph && this->index()==n.index()){
        return true;
      }else{
        return false;
      }
      //(void) n;          // Quiet compiler warning
      //return false;
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
      //If one node is in a smaller graph than the other or 
      //index of one node is smaller than the other
      if ((this->graph<n.graph)||(this->index()<n.index())){
        return true;
      }else{
        return false;
      }
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    //Pointer back tp Graph container
    Graph* graph;
    //index of this node
    size_type node_idx;
    //Constructor to construct Node objects from Graph class methods
    Node(const Graph* new_graph, size_type new_node_idx)
        : graph(const_cast<Graph*>(new_graph)),node_idx(new_node_idx){
        }

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  std::map<size_type,std::tuple<Point*,node_value_type>> get_nodes()const{
    return this->nodes;
  }
  std::map<size_type,std::vector<size_type>> get_edge_list()const{
    return this->edge_list_;
  }
  size_type size() const {
    return this->nodes_.size();
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
    // HW0: YOUR CODE HERE
    //dynamically allocate position
    Point* pos = new Point;
    *pos = position;
    size_type new_node_idx = this->size();
    std::tuple<Point*,node_value_type> new_node(pos,node_value_type());
    this->nodes_[new_node_idx] = new_node;
    //check if num_nodes was modified correctly
    assert(this->num_nodes() == new_node_idx +1);
    std::vector<size_type> new_vec;
    this->edge_list_[new_node_idx] = new_vec;
    return Node(this,new_node_idx);
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph==this && n.index()<this->size()){
      return true;
    }else{return false;}
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
    assert(0<=i && i<num_nodes());
    return Node(this,i);

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
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      //return Node();      // Invalid Node
      return Node(this->graph,node1_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph,node2_idx);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
        //If first node are the same and second node are the same.
        //or they are the same but in different order
        //Even though we don't add the edge if the same edge but different
        //node order exists, maybe we should still make this method able
        //to compare for this situation...
        if(((this->node1() == e.node1()) && (this-> node2() == e.node2()))
        ||((this->node1() == e.node2()) && (this-> node2() == e.node1()))){
        return true;
      }
      else{return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    std::vector<Node> node_ordered(const Edge& edge0) const {
        std::vector<Node> ordered_nodes;
        if(edge0.node1()<edge0.node2()){
          ordered_nodes.push_back(edge0.node1());
          ordered_nodes.push_back(edge0.node2());
        }else{ordered_nodes.push_back(edge0.node2());
          ordered_nodes.push_back(edge0.node1());}
        return ordered_nodes;
      }

    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //Same as operator==, we consider ordering in this method
        std::vector<Node> order1 = node_ordered(*this);
        std::vector<Node> order2 = node_ordered(e);

        if((order1[0]<order2[0])||((order1[0]==order2[0])&&(order1[1]<order2[1]))){
              return true;
          }else{return false;}
      }
      


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type node1_idx;
    size_type node2_vec_idx;
    size_type node2_idx;
    Edge(const Graph* new_graph,size_type new_node1_idx,size_type new_node2_vec_idx):
      graph(const_cast<Graph*>(new_graph)),node1_idx(new_node1_idx),node2_vec_idx(new_node2_vec_idx){
        if(node1_idx<graph->edge_list_.size()){
        node2_idx = graph->edge_list_.at(node1_idx)[node2_vec_idx];
        }
        else{node2_idx = 0;}
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type num_edges = 0;
    for (auto i: edge_list_){
        num_edges = num_edges + i.second.size();
    }
    num_edges = num_edges/2;
    return num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
    assert(0<=i && i<this->num_edges());
    //use iterator to return the ith edge
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;
    assert(a.graph==this);
    assert(b.graph==this);
    size_type a_idx = a.index();
    size_type b_idx = b.index();
    //check in both orders
    std::vector<size_type> a_edges = edge_list_.at(a_idx);
    for (auto v: a_edges){
      if(v==b_idx){
        return true;
      }else{continue;}
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
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
    // assert(a.graph==this && b.graph==this);
    assert(a.index() != b.index());
    // assert(a.index()<this->size() && b.index()<this->size());
    // size_type old_num_edges = num_edges();
    std::vector<size_type> a_vec = edge_list_.at(a.index());

    for (size_type i = 0;i<a.degree();++i){
      if (a_vec[i] == b.index()){
        Edge new_edge(this,a.index(),i);
        // assert taking too long time...
        // assert(old_num_edges == num_edges());
        return new_edge;
      }
    }
    //add new edge to edge_list_
    this->edge_list_[a.index()].push_back(b.index());
    this->edge_list_[b.index()].push_back(a.index());
    // check if correctly added
    assert(has_edge(a,b));
    assert(has_edge(b,a));
    //check new num_edges
    // assert(old_num_edges+1 == num_edges());
    size_type b_in_vec = edge_list_[a.index()].size()-1;
    return Edge(this,a.index(),b_in_vec);
    

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    
    //delete pointers
    for (auto p: this->nodes_){
      delete p;
    }
    this->nodes_.clear();
    this->edge_list_.clear();
    assert(num_nodes()==0);
    assert(num_edges()==0);

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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator():graph(nullptr) {
    }
  
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    Node operator*() const{
      return Node(graph,index);
    }
    NodeIterator& operator++(){
      //When index too large or at the end, set it to size of nodes_
      if (index >= graph_size-1){
        index = graph_size;
        return *this;
      }else{
        index = index+1;
        return *this;
      }
    
    }
    bool operator==(const NodeIterator& iter) const{
      return ((this->graph == iter.graph) && (this->index == iter.index));
    }
    
    private:
    friend class Graph;
    const Graph* graph;
    size_type index;
    size_type graph_size;
    NodeIterator(const Graph* g, size_type idx):graph(g),index(idx){graph_size = graph->num_nodes();}
    // HW1 #2: YOUR CODE HERE
    
  };
  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }
  node_iterator node_end() const{
    return NodeIterator(this,this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator():graph(nullptr) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
          //return the edge, constructor pass in vector index of node2, not node index
           return Edge(graph,node1_idx,node2_idx);
    }
    IncidentIterator& operator++(){
      //If given node2_idx at the end of vector or probably too large
      //set it to vector size of edge_list_[node1_idx] to prevent seg fault
      size_type vec1_size = graph->edge_list_.at(node1_idx).size();
      if (node2_idx >= vec1_size-1){
        node2_idx = vec1_size;
        return *this;
      }else{
        node2_idx = node2_idx+1;
        return *this;
      }
    }
    bool operator==(const IncidentIterator& iter) const{
      return ((iter.graph == this->graph)&&(iter.node1_idx == this->node1_idx) 
                  && (iter.node2_idx == this->node2_idx));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph;
    size_type node1_idx;
    size_type node2_idx;
    IncidentIterator(const Graph* g,size_type node1, size_type node2):
                      graph(g),node1_idx(node1),node2_idx(node2){}
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
    EdgeIterator():graph(nullptr){
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      //Return the edge, notice Edge was constructed by vector index of node2
      //instead of node index
      return Edge(graph,node1_idx,node2_idx);
    }
    EdgeIterator& operator++(){
        //When node2_idx is not at the end of its vector, first increment inside this vector
        //When at the end, move to next node1_idx
        while(node1_idx<vec1_size){
          vec2_size = graph->edge_list_.at(node1_idx).size();
          while(node2_idx<vec2_size-1){
            node2_idx = node2_idx+1;
            //To only visit each edge once, only return the iterator when node1_idx is
            //smaller than the node index of node2
            if (node1_idx < graph->edge_list_.at(node1_idx)[node2_idx]){
              return *this;
            }
          }
          node1_idx = node1_idx+1;
          node2_idx = 0;
        }
        //If haven't return after while loop, this means the iterator has reach the end
        node1_idx = vec1_size;
        node2_idx = 0;
        return *this;

        
      }
      

    
    bool operator==(const EdgeIterator& iter) const{
      return ((this-> graph == iter.graph)&&(this->node1_idx== iter.node1_idx)
              &&(this->node2_idx== iter.node2_idx));
    }

   private:
    friend class Graph;
    const Graph* graph;
    size_type node1_idx;
    size_type node2_idx;
    size_type vec1_size;
    size_type vec2_size;
    //EdgeIterator Constructor: when passing in node1_idx too large, set vec2_size to 0 to prevent seg fault
    EdgeIterator(const Graph* g,size_type node1,size_type node2):graph(g),node1_idx(node1),node2_idx(node2){
                          vec1_size= graph->edge_list_.size();
                          if (node1_idx >= vec1_size){
                            vec2_size = 0;
                          }else{
                          vec2_size= graph->edge_list_.at(node1_idx).size();};
                        }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  public:
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0,0);
  }
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->edge_list_.size(),0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
