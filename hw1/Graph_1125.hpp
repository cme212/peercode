#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <string>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

  // ===========================================================================
  // GRAPH
  // ===========================================================================

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 *
 * This code uses the proxy design pattern: the Node and Edge classes are 
 * used as proxies, while the internal_node and internal_edge structs contain
 * the underlyind data of those structures, allowing better data encapsulation.
 * 
 * The code has two functions, fetch_node() and fetch_edge(), that allow to 
 * access the data in the internal structs when we have a Node or Edge object.
 */

template <typename V>
class Graph: private totally_ordered <Graph<V>> {

  // Predeclare the internal struct for node and eges
  struct internal_edge; 
  struct internal_node;

public:

  // ===========================================================================
  // PUBLIC TYPE DEFINITIONS
  // ===========================================================================

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

  using node_value_type = V;

  // ===========================================================================
  // CONSTRUCTORS AND DESTRUCTOR
  // ===========================================================================

  /** Constructs an empty graph. */
  Graph() : nodes_(), edges_(), internal_nodes_(), internal_edges_(), 
  nodes_to_edge_() {};

  /** Default destructor */
  ~Graph() = default;

  // ===========================================================================
  // NODES
  // ===========================================================================

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered <Node> {
   
   public:
    /** Constructs an invalid node. */
    Node() {}

    /** Returns this node's position (cannot be modified). */
    const Point& position() const 
    {
      return *(fetch_node().pt_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const 
    {
      return node_id_;
    }

    IncidentIterator incident_begin()
    {
      return IncidentIterator(this->graph_, get_incident_edges(), 0);
    }

    IncidentIterator incident_end()
    {
      return IncidentIterator(this->graph_, get_incident_edges(), degree());
    }

    /** Return the number of incident edges of a node. **/
    size_type degree() const 
    {
      std::vector<edge_type> incident_edges = get_incident_edges();
      return incident_edges.size();
    }

    /** Return a vector containing the incident edges of this node **/
    std::vector<edge_type> get_incident_edges() const 
    {
      std::vector<edge_type> incident_edges;
      EdgeIterator e_it = (*(this->graph_)).edge_begin();

      while(e_it != (*(this->graph_)).edge_end())
      {
        Edge e = *e_it;
        Node n1 = e.node1();
        Node n2 = e.node2();
        if(n1 == (*this) || n2 == (*this))
          incident_edges.push_back(e);
        ++e_it;
    }
      return incident_edges;
    }

    /** Return this node's value */
    node_value_type value()
    {
      return fetch_node().get_node_value();
    }

    /** Return this node's value (cannot be modified) */
    const node_value_type value() const
    {
      return fetch_node().get_node_value();
    }

    /** Set the value of an internal node. */
    void set_node_val(node_value_type val) 
    {
      fetch_node().set_node_value(val);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const 
    {
      if(n.node_id_ == node_id_ && n.graph_ == graph_)
        return true;
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
    bool operator<(const Node& n) const 
    {
      if(n.node_id_ > node_id_ && n.graph_ == graph_)
        return true;
      return false;
    }

    /** Helper method to return the appropriate element.
     * This loops over the elements until it finds the element with the
     * correct uid.
     */
    internal_node& fetch_node() const 
    {
        return this->graph_->internal_nodes_.find(node_id_)->second;
    }

    private:
      friend class Graph;
      Graph* graph_; // Pointer to the graph of the node
      size_type node_id_; // Index of the node

      /** Construct a valid node using initalizer list. */
      Node(const Graph* graph, size_type idx):
        graph_(const_cast<Graph*>(graph)), node_id_(idx){}
  };

  // ===========================================================================
  // GRAPH - NODES METHOD
  // ===========================================================================

  void default_initialize_nodes(int a) const
  {
    // Declare the iterators
    node_iterator niter_first = node_begin();

    // Default initialize values of all nodes to -1.
    for(; niter_first != node_end(); ++niter_first)
      (*niter_first).set_node_val(a);
  }

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const 
  {
      return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const 
  {
    return size();
  }

  Node add_node(const Point& position, const node_value_type& a = node_value_type()) 
  {
    // Declare a new point object and make it point to the input address
    Point* new_point = new Point;
    *new_point = position;

    // Declare a new node and add it to the vector of nodes, Create the pair for the map
    internal_node new_internal_node = internal_node(new_point, a); 

    // Create the pair for the map
    std::pair<size_type, internal_node> key_val (num_nodes(), new_internal_node);
    internal_nodes_.insert(key_val);

    // Declare new node and push it to vector
    node_type new_node = Node(this, num_nodes()); 
    nodes_.push_back(new_node);

    // Return the new node using the current graph and index
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const 
  {  
    if(this == n.graph_ && n.node_id_ < size())
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const 
  {
    return Node(this, i);
  }

  std::string create_node_id(const Node&a, const Node&b) const
  {
    // Get the nodes ids
    std::string a_id = std::to_string(a.node_id_);
    std::string b_id = std::to_string(b.node_id_);

    // Create key
    std::string key = a_id + "-" + b_id;

    return key;
  }

  // ===========================================================================
  // EDGES
  // ===========================================================================

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:

    /** Construct an invalid Edge. */
    Edge() {}

    /** Returns the first node of this Edge */
    Node node1() const 
    {
      return Node(this->graph_, fetch_edge().first_node_idx_);      
    }

    /** Return the second node of this Edge */
    Node node2() const 
    {
      return Node(this->graph_, fetch_edge().second_node_idx_);      
    }

    /** Returns the first node of this Edge */
    size_type node1_idx() const 
    {
      return fetch_edge().first_node_idx_;      
    }

    /** Returns the first node of this Edge */
    size_type node2_idx() const 
    {
      return fetch_edge().second_node_idx_;      
    }

    size_type get_edge_id() const 
    {
      return edge_id_;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * First, check if both edges are from same graph. Then, check if 
     * their nodes match (in one direction or the other as graph is undirected)
     */
    bool operator==(const Edge& e) const 
    {
      if((this->graph_ == e.graph_) && ((fetch_edge().first_node_idx_== 
        e.fetch_edge().first_node_idx_ && fetch_edge().second_node_idx_ == 
        e.fetch_edge().second_node_idx_) || (fetch_edge().first_node_idx_ == 
        e.fetch_edge().second_node_idx_ && fetch_edge().second_node_idx_
        == e.fetch_edge().second_node_idx_)))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const 
    {
      assert(e.graph_ != NULL and this->graph_ != NULL);

      // Get the min and max of nodes indexes for the two edges 
      size_type e_min = std::min(e.fetch_edge().first_node_idx_, 
        e.fetch_edge().second_node_idx_);
      size_type e_max = std::max(e.fetch_edge().first_node_idx_, 
        e.fetch_edge().second_node_idx_);

      size_type this_min = std::min(fetch_edge().first_node_idx_, 
        fetch_edge().second_node_idx_);
      size_type this_max = std::max(fetch_edge().first_node_idx_, 
        fetch_edge().second_node_idx_);

      if ((graph_ < e.graph_) ||((graph_ == e.graph_) && ((this_min < e_min) ||
       ((this_min == e_min) && (this_max < e_max))))) 
      {
        return true;
    }
    return false;
    }

    // Fetching function to get the underlying data of an Edge
    internal_edge& fetch_edge() const 
    {
      return this->graph_->internal_edges_.find(edge_id_)->second;
    }

    private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      // Edge's private data attributes
      Graph* graph_;
      size_type edge_id_;

      //Construct a valid Edge using initalizer list.
      Edge(const Graph* graph, size_type edge_id):
        graph_(const_cast<Graph*>(graph)),
        edge_id_(edge_id){};
  };

  // ===========================================================================
  // GRAPH - EDGES METHOD
  // ===========================================================================

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const 
  {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const 
  {
    return Edge(this, i);
  }

  Edge find_edge(const Node& a, const Node& b) 
  {
    std::string key1 = create_node_id(a, b);
    std::string key2 = create_node_id(b, a);

    if(nodes_to_edge_.find(key1) == nodes_to_edge_.end()) 
    {
      size_type edge_num = nodes_to_edge_.find(key2)->second;
      return edge(edge_num);
    }
    else if(nodes_to_edge_.find(key2) == nodes_to_edge_.end()) 
    {
      size_type edge_num = nodes_to_edge_.find(key1)->second;
      return edge(edge_num);
    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) 
  {
    std::string key1 = create_node_id(a, b);
    std::string key2 = create_node_id(b, a);

    if(nodes_to_edge_.find(key1) == nodes_to_edge_.end() &&
      nodes_to_edge_.find(key2) == nodes_to_edge_.end())
        return false;
    else
      return true;
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
  Edge add_edge(const Node& a, const Node& b) 
  {
    // Check if the edge already exists. If yes, return it
    if(has_edge(a, b))
      return find_edge(a, b);

    // Create a key, val for the unordered map
    std::string key = create_node_id(a, b);
    std::pair<std::string, size_type> key_val (key, num_edges());
    nodes_to_edge_.insert(key_val);

    // Create a new edge to be added
    internal_edge new_internal_edge = internal_edge(a.node_id_, b.node_id_);

    // Create the pair for the map
    std::pair<size_type, internal_edge> key_val_2 (num_edges(), new_internal_edge);
    internal_edges_.insert(key_val_2);

    // Declare new edge and push it to vector
    edge_type new_edge = Edge(this, num_edges()); 
    edges_.push_back(new_edge);

    // Return the new Edge
    return Edge(this, num_edges() - 1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() 
  {
    nodes_.clear(); edges_.clear();
    internal_nodes_.clear();internal_edges_.clear(); nodes_to_edge_.clear();
  }

  // ===========================================================================
  // NODE ITERATOR
  // ===========================================================================

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
    NodeIterator(){};

    NodeIterator operator++()
    {
      node_iter_idx_++;
      return (*this);
    }

    bool operator==(const NodeIterator& n) const
    {
      if(this->node_iter_idx_ == n.node_iter_idx_ 
        && this->graph_iter_ == n.graph_iter_)
        return true;
      return false;
    }

    //Defines inequality between two iterators
    bool operator!=(const NodeIterator& n) const
    {
        // Reuse the == operator
        return !(n == (*this));
    }

    Node operator*() const
    {
      return Node(graph_iter_, node_iter_idx_);
    }

    size_type get_node_iter_idx() const
    {
      return node_iter_idx_;
    }

    NodeIterator(const Graph* graph, size_type node_iter_idx):
    node_iter_idx_(node_iter_idx),
    graph_iter_(const_cast<Graph*>(graph)) {};

   private:
    friend class Graph;
    size_type node_iter_idx_;
    Graph* graph_iter_; // Pointer to the graph of the node
  };

  // ===========================================================================
  // GRAPH - NODE ITERATOR METHODS
  // ===========================================================================

    NodeIterator node_begin() const
    {
      return NodeIterator(const_cast<Graph*>(this), 0);
    }

    NodeIterator node_end() const
    {
      return NodeIterator(const_cast<Graph*>(this), nodes_.size());
    }

  // ===========================================================================
  // INCIDENT ITERATOR
  // ===========================================================================

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
    IncidentIterator() {};

    Edge operator*() const 
    {
      return Edge(graph_ptr_inc_itr_, incident_edges_.at(edge_iter_incident_idx_).get_edge_id());
    }

    IncidentIterator& operator++() 
    {
      edge_iter_incident_idx_++;
      return (*this);
    }

    bool operator==(const IncidentIterator& i) const 
    {
      if(this->graph_ptr_inc_itr_ == i.graph_ptr_inc_itr_ 
        && this->incident_edges_ == i.incident_edges_
        && this->edge_iter_incident_idx_ == i.edge_iter_incident_idx_)
        return true;
      return false;
    }

    //Defines inequality between two iterators
    bool operator!=(const IncidentIterator& i) const
    {
        // Reuse the == operator
        return !(i == (*this));
    }

    // Get incident edges vector of an IncidentIterator
    std::vector<edge_type> get_incident_edges_vec() 
    {
      return incident_edges_;
    }

   private:
    friend class Graph;
    friend class Node;
    Graph* graph_ptr_inc_itr_;
    std::vector<edge_type> incident_edges_;
    size_type edge_iter_incident_idx_;

    IncidentIterator(
      const Graph* graph, 
      std::vector<edge_type> incident_edges,
      size_type edge_iter_incident_idx):
      graph_ptr_inc_itr_(const_cast<Graph*>(graph)),
      incident_edges_(incident_edges),
      edge_iter_incident_idx_(edge_iter_incident_idx)
      {};
  };

  // ===========================================================================
  // EDGE ITERATOR
  // ===========================================================================

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
    EdgeIterator() {}

    EdgeIterator& operator++() 
    {
      edge_iter_idx_++;
      return (*this);
    }

    bool operator==(const EdgeIterator& n) const
    {
      if(this->edge_iter_idx_ == n.edge_iter_idx_ 
        && this->graph_iter_ == n.graph_iter_)
        return true;
      return false;
    }

    bool operator!=(const EdgeIterator& e) const
    {
        // Reuse the == operator
        return !(e == (*this));
    }

    Edge operator*() const
    {
      return Edge(graph_iter_, edge_iter_idx_);
    }

   private:
    friend class Graph;
    size_type edge_iter_idx_;
    Graph* graph_iter_; // Pointer to the graph of the edge

    EdgeIterator(Graph* graph, size_type edge_iter_idx):
    edge_iter_idx_(edge_iter_idx),
    graph_iter_(graph) {};

  };

  // ===========================================================================
  // GRAPH - EDGE ITERATOR METHODS
  // ===========================================================================

    EdgeIterator edge_begin()
    {
      return EdgeIterator(this, 0);
    }

    EdgeIterator edge_end()
    {
      return EdgeIterator(this, edges_.size());
    }

private:

  // ===========================================================================
  // INTERNAL NODE
  // ===========================================================================
  
    struct internal_node 
    {
    // Data attributes of the internal_node struct: a Point object and unique id
    const Point* pt_;
    node_value_type value_;

    /** Change this node's value */
    void set_node_value(node_value_type val)
    {
      value_ = val;
    }

    /** Return this node's value */
    node_value_type get_node_value()
    {
      return value_;
    }

    // Default constructor
    internal_node(): pt_(nullptr), value_(-1){};

    // Parameterized Constructor (with value)
    internal_node(const Point* position, node_value_type value): 
    pt_(const_cast<Point*>(position)), 
    value_(value) {};
  };

  // ===========================================================================
  // INTERNAL EDGE
  // ===========================================================================

    struct internal_edge 
    {
    // Data attributes of the internal_edge struct: a unique edge idx
    // the ids of the nodes it is connecting
    size_type first_node_idx_;
    size_type second_node_idx_;

    /** Return this edge's first node index */
    size_type get_edge_n1()
    {
      return first_node_idx_;
    }

    /** Return this edge's second node index */
    size_type get_edge_n2()
    {
      return second_node_idx_;
    }
    // Default constructor
    internal_edge() {};

    // Parameterized Constructor
    internal_edge(size_type first_idx, size_type second_idx): 
    first_node_idx_(first_idx), 
    second_node_idx_(second_idx)
    {};
  };

  // ===========================================================================
  // GRAPH DATA ATTRIBUTES
  // ===========================================================================

  // Vectors to store the Nodes and Edges in the Graph
  std::vector<node_type> nodes_;
  std::vector<edge_type> edges_;

  // Maps to store the Internal Nodes and Edges in the Graph
  std::map<size_type, internal_node> internal_nodes_;
  std::map<size_type, internal_edge> internal_edges_;

  // Unordered map to store the Internal Nodes and Edges in the Graph
  std::unordered_map<std::string, size_type> nodes_to_edge_;
};

#endif // CME212_GRAPH_HPP
