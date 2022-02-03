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
template <typename V>
class Graph
{
private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  unsigned node_size_;
  unsigned edge_size_;
  unsigned edge_counts;
  struct proxyNode;
  struct proxyEdge;
  std::vector<proxyNode> Nodes;
  std::vector<proxyEdge> Edges;
  std::map<std::vector<unsigned>, unsigned> e2ID;

public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
      

  /** Predeclaration of Node type. */
  //class Node;
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
  Graph()
  {
    node_size_ = 0;
    edge_size_ = 0;
    edge_counts=0;
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
  class Node : private totally_ordered<Node>
  {
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
      //this->graph_ = nullptr;
      //this->id=0;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point &position() const
    {
      // HW0: YOUR CODE HERE
      return graph_->Nodes[id].pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return this->id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
     //size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    node_value_type &value(){
      return graph_->Nodes[id].val;
    }
    const node_value_type &value() const{
      return graph_->Nodes[id].val;
    }
    size_type degree() const{
      return graph_->Nodes[id].neighbors.size();
    }
    incident_iterator edge_begin() const{
      return incident_iterator(this->graph_,this->id,0);
    }
    incident_iterator edge_end() const{
      return incident_iterator(this->graph_, this->id, this->degree());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      // HW0: YOUR CODE HERE
      //(void)n; // Quiet compiler warning
      return this->graph_==n.graph_ and this->id == n.id;
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
    bool operator<(const Node &n) const
    {
      // HW0: YOUR CODE HERE
     // Quiet compiler warning
      return this->id < n.id;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph *graph_;
    size_type id;
    Node( const Graph *graph, size_type id): graph_{const_cast<Graph*>(graph)},id{id}{}
   
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    // HW0: YOUR CODE HERE
    return node_size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type & v = node_value_type())
  {
    // HW0: YOUR CODE HERE
    std::vector<size_type> neighbors=std::vector<size_type>();
    node_size_ += 1;
    proxyNode proxynode = proxyNode(node_size_ - 1, position,v,neighbors);
    Nodes.push_back(proxynode);
    return Node(this, node_size_ - 1);

    // return Node(); // Invalid node
  }
  //Node add_node(const Point &, const node_value_type & = node_value_type());

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    // HW0: YOUR CODE HERE
    //(void)n; // Quiet compiler warning
    return n.id < node_size_ and n.id >= 0;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    // HW0: YOUR CODE HERE
    //(void)i;       // Quiet compiler warning
    if (i >= node_size_)
    {
      return Node();
    }
    else
    {
      return Node(this, i);
    }
    // Invalid node
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
  class Edge : private totally_ordered<Node>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      // HW0: YOUR CODE HERE
      size_type node_id = graph_->Edges[this->id].node1;

      return graph_->node(node_id); // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      size_type node_id = graph_->Edges[this->id].node2;
      return graph_->node(node_id); // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      (void)e; // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return this->id == e.id;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      (void)e; // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return this->id < e.id;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph *graph_;
    size_type id;
    Edge( const Graph *graph, size_type id): graph_{const_cast<Graph* const>(graph)},id{id}{}
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    // HW0: YOUR CODE HERE
    return this->edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    // HW0: YOUR CODE HERE
    (void)i;    
    return Edge(this, i); // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    // HW0: YOUR CODE HERE
    if (!has_node(a) or !has_node(b))
    {
      return false;
    }
    std::vector<size_type> node_set = {a.id, b.id};

    auto search = e2ID.find(node_set);
    //std::vector<size_type> node_set2 = {b.id, a.id};
    //auto search2 = e2ID.find(node_set2);

    if (search != e2ID.end())
    {
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
  Edge add_edge(const Node &a, const Node &b)
  {
    std::vector<size_type> key = {a.id, b.id};
    std::vector<size_type> key2 = {b.id, a.id};

    if (has_edge(a, b))
    {
      size_type id = e2ID.find(key)->second;
      return Edge(this, id);
    }
    //(void)a, (void)b; // Quiet compiler warning
    else{
    this->edge_size_ += 1;
    this->edge_counts+=2;

    proxyEdge proxyedge = proxyEdge(a.id, b.id, edge_counts - 2);
    e2ID[key] = edge_counts - 2;

    proxyEdge proxyEdge2 = proxyEdge(b.id, a.id, edge_counts - 1);
    e2ID[key2] = edge_counts - 1;

    Edges.push_back(proxyedge);
    Edges.push_back(proxyEdge2);
   

    Nodes[a.id].neighbors.push_back(b.id);
    Nodes[b.id].neighbors.push_back(a.id);
    
    return Edge(this, edge_counts - 2); 
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    Nodes.clear();
    Edges.clear();
    e2ID.clear();
    node_size_ = 0;
    edge_size_ = 0;
    edge_counts=0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::forward_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }


    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Node operator*() const{

       return this->graph_->node(this->id);

     }
      NodeIterator& operator++(){
       if (this->id == this->graph_->size()){
         return *this;
       }
       else{
         this->id++;
         return *this;
     }
      }
     bool operator==(const NodeIterator& node) const{
       return this->id == node.id;
     }

     bool operator!=(const NodeIterator &node) const
     {
      return this->id != node.id;
     }

  private:
    friend class Graph;
    Graph *graph_;
    size_type id;

    NodeIterator( const Graph *graph, size_type id_) : graph_{const_cast<Graph*>(graph)} ,id{id_} {}
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const
  {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const
  {
    return NodeIterator(this, this->size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
      this->graph_=nullptr;
      this->id=0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
    Node a= Node(this->graph_,this->Node_id);
    Node b= Node(this->graph_,graph_->Nodes[this->Node_id].neighbors[this->id]);
    std::vector<unsigned int> key = {a.id,b.id};
  
    size_type edge_id=this->graph_->e2ID.find(key)->second;
    return Edge(this->graph_,edge_id);
    }

    IncidentIterator& operator++(){
     
      
        this->id++;
        return *this;
      
    }
    bool operator==(const IncidentIterator& inc) const{
      return this->Node_id == inc.Node_id and this->id==inc.id;
    }

    bool operator!=(const IncidentIterator &inc) const
    {
      return this->Node_id != inc.Node_id or this->id != inc.id;
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
     Graph* graph_;
     size_type Node_id;
     size_type id;
     IncidentIterator(const Graph* graph,int Node_id_,size_type id_):graph_(const_cast<Graph*>(graph)),Node_id(Node_id_),id(id_){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {

    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Edge operator*() const{
       size_type edge_id= this->graph_->Edges[this->id].index;
       return Edge(this->graph_,edge_id);
     }
     EdgeIterator& operator++(){
       //if (this->id<this->graph_->size()){
         //return *this;
       //}
       //else{
        
         this->id+=2;
         //this->id++;
         //return *this;
     //  }
        return *this;
     }
     bool operator==(const EdgeIterator& edgeit) const{
       return this->id == edgeit.id and this->graph_==edgeit.graph_;
     }

     bool operator!=(const EdgeIterator &edgeit) const
     {
       return this->id != edgeit.id or this->graph_ != edgeit.graph_;
     }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type id;
    EdgeIterator(const Graph* graph, size_type id_):graph_{const_cast<Graph*>(graph)},id{id_}{}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const{
     return EdgeIterator(this,0);
   }
   edge_iterator edge_end() const{
     return EdgeIterator(this,this->edge_counts);
   }

private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct proxyNode
  {
    size_type index;
    Point pt;
    node_value_type val;
    std::vector<size_type> neighbors;
    proxyNode(size_type index, Point pt, node_value_type val,std::vector<size_type> neighbors)
    {
      this->index = index;
      this->pt = pt;
      this->val=val;
      this->neighbors=neighbors;
    }
  };

  struct proxyEdge
  {
    size_type index;
    size_type node1;
    size_type node2;
    proxyEdge(size_type node1, size_type node2, size_type index)
    {
      this->node1 = node1;
      this->node2 = node2;
      this->index = index;
    }
  };
};

#endif // CME212_GRAPH_HPP
