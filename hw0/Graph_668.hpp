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
 *
 * This code uses the proxy design pattern: the Node and Edge classes are 
 * used as proxies, while the internal_node and internal_edge structs contain
 * the underlyind data of those structures, allowing better data encapsulation.
 * 
 * The code has two functions, fetch_node() and fetch_edge(), that allow to 
 * access the data in the internal structs when we have a Node or Edge object.
 */

class Graph {

  // Predeclare the internal struct for node and eges
  struct internal_edge; 
  struct internal_node;

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
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Constructs an empty graph. */
  Graph() : nodes_(std::vector<internal_node>(0)), 
    edges_(std::vector<internal_edge>(0)){}

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
  class Node {
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
      return fetch_node().node_idx_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const 
    {
      if(n.node_id_ == fetch_node().node_idx_ && n.graph_ == graph_)
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
      if(n.node_id_ > fetch_node().node_idx_ && n.graph_ == graph_)
        return true;
      return false;
    }

    private:
      friend class Graph;
      Graph* graph_; // Pointer to the graph of the node
      size_type node_id_; // Index of the node

      /** Construct a valid node using initalizer list. */
      Node(const Graph* graph, size_type idx):
        graph_(const_cast<Graph*>(graph)), node_id_(idx){}

      /** Helper method to return the appropriate element.
       * This loops over the elements until it finds the element with the
       * correct uid.
       */
      internal_node fetch_node() const {
        for (size_type i = 0; i < graph_->nodes_.size(); ++i)
          if (graph_->nodes_[i].node_idx_ == node_id_)
            return graph_->nodes_[i];
        assert(false);
      }
  };

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

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) 
  {
    // Declare a new point object and make it point to the input address
    Point* new_point = new Point;
    *new_point = position;

    // Declare a new node and add it to the vector of nodes
    internal_node new_node = internal_node(new_point, size()); 
    nodes_.push_back(new_node);

    // Return the new node using the current graph and index
    return Node(this,size()-1);
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

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
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
       ((this_min == e_min) && (this_max < e_max))))) {
        return true;
    }
    return false;
    }

    private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      // Edge's private data attributes
      Graph* graph_;
      size_type edge_id_;

      //Construct a valid Edge using initalizer list.
      Edge(const Graph* graph,size_type edge_id):
        graph_(const_cast<Graph*>(graph)),edge_id_(edge_id){}
      
      // Fetching function to get the underlying data of an Edge
      internal_edge fetch_edge() const 
      {
        for (size_type i = 0; i < graph_->edges_.size(); ++i)
          if (graph_->edges_[i].edge_idx_ == edge_id_)
            return graph_->edges_[i];
        assert(false);
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const 
  {
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const 
  {
    return Edge(this,i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const 
  {
    assert(this == a.graph_ && this == b.graph_);

    for(size_type i = 0;i<this->edges_.size();++i)
      if((edges_[i].first_node_idx_ == a.node_id_ && 
        edges_[i].second_node_idx_ == b.node_id_) || 
        (edges_[i].first_node_idx_ == b.node_id_ &&
         edges_[i].second_node_idx_ == a.node_id_))
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

  Edge add_edge(const Node& a, const Node& b) 
  {
    // Check if the edge already exists. If yes, return it
    if(has_edge(a,b))
      return Edge(this, num_edges());

    // If not, create a new edge to be added
    internal_edge new_edge = internal_edge(a.node_id_, b.node_id_, num_edges());
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

   private:
    friend class Graph;
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

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

private:

    struct internal_edge 
    {
    // Data attributes of the internal_edge struct: a unique edge idx
    // the ids of the nodes it is connecting
    size_type first_node_idx_;
    size_type second_node_idx_;
    size_type edge_idx_;

    // Default constructor
    internal_edge() {};

    // Parameterized Constructor
    internal_edge(size_type first_idx, size_type second_idx, size_type edge_idx)
    : first_node_idx_(first_idx), second_node_idx_(second_idx), 
    edge_idx_(edge_idx){};

  };
    struct internal_node 
    {
    // Data attributes of the internal_node struct: a Point object and unique id
    const Point* pt_;
    size_type node_idx_;

    // Default constructor
    internal_node(): pt_(nullptr), node_idx_(-1) {};

    // Parameterized Constructor
    internal_node(const Point* position, size_type id): pt_(const_cast<Point*>(position)),
    node_idx_(id){};

    };

  //Vectors to store the Nodes and Edges in the Graph
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

};

#endif // CME212_GRAPH_HPP
