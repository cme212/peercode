#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
  //class internal_node; // Predeclaring the internal struct 
  class internal_edge;
  //struct HashEdge; 
  //struct NodeSet;


 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


 public:
  class internal_node;

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node; // node_type is an alias for Node

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
  using size_type = unsigned; //Short for unsigned int
  using key_pair_type = std::tuple<size_type, size_type>;
  
  std::unordered_map<size_type, size_type> uid_to_idx;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() // A way of constructing an empty graph 
    : nodes_(), edges_() { //Technically we don't need to populate it
  }
    /** TODO: Consider adding edges_() to the above, if I decide to make it a proxy class */

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   * The Node class is a proxy class for the points
   */
  class Node {
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
    Node() { // Creates invalid simple element
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {// HW0: YOUR CODE HERE
      //TODO: Consider adding a map between uids and idx   
      size_type get_idx = graph_->uid_to_idx[uid_]; //Get the idx
      return graph_->nodes_[get_idx].point; //
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { // HW0: YOUR CODE HERE
      //size_type get_idx = graph_->; //Get the idx
      return uid_;
      //return graph_->nodes_[get_idx].uid;
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
    bool operator==(const Node& n) const { // HW0: YOUR CODE HERE
      // If nodes has the same id, then true. Else, false.
     if (graph_->nodes_[uid_].uid == n.uid_) {
       return true;
     } else {
       return false;
     }
    //  (void) n;          // Quiet compiler warning
    //  return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const { // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].uid < n.uid_; 
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private: //Private attributes of node
    // Pointer back to the graph containeer;
    Graph* graph_; 
    // This node's unique identification number
    size_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    friend class Graph;
    friend class internal_node;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  }; //This is the end of the proxy class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { // HW0: YOUR CODE HERE
    return nodes_.size(); // Return the size attribute
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size(); //It calls the size() method
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point position) {// HW0: YOUR CODE HERE
    uid_to_idx[num_nodes()] = num_nodes(); // 
    nodes_.push_back(internal_node(position, num_nodes())); //Add a new node instance to the vector
    //std::cout << "This node's position is: " << position << std::endl;
    // Returns a Node that points to the new node
    return Node(this, num_nodes()-1);        // Invalid node
    //
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const { // HW0: YOUR CODE HERE
    //std::cout << "The n.index is: " << n.index() << " and the num_nodes is: " << num_nodes() << std::endl;
    if ((n.index() < num_nodes()) && (n.index() >= 0)) {
      return true;
    } else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { // HW0: YOUR CODE HERE
    //std::cout << "number of nodes are:" << num_nodes() << "and the index we want is: " << i << std::endl;
    assert(i < num_nodes()); // Check that the idx is within the size of the graph
    assert(i >= 0);
    return Node(this, i); // Return the value of the internal node ID
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
    Edge() { // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const { // HW0: YOUR CODE HERE
      return graph_->node(a_idx);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const { // HW0: YOUR CODE HERE
      return graph_->node(b_idx);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      if ((e.node1().index() == a_idx && e.node2().index() == b_idx) || (e.node1().index() == b_idx && e.node2().index() == a_idx)) {
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const { //HW0: YOUR CODE HERE
      if (eid_ < e.eid_) {
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Creating a graph pointer
    Graph* graph_; 

    size_type eid_, a_idx, b_idx; //The two nodes connected to the edge

    /** Private Constructor */
    Edge(size_type eid, size_type a_idx, size_type b_idx) 
      : eid_(eid), a_idx(a_idx), b_idx(b_idx) { //Each edge will have node a_ and node b_
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { // HW0: YOUR CODE HERE
    return edges_.size(); //Return the number of edges in the container
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { // HW0: YOUR CODE HERE
    assert(i >= 0); //Check that i is greater than or equal to zero
    assert(i < num_edges()); // Checkt that i is less than the number of edges in the vector
    //Need the edge id, a, and b
    size_type id = edges_[i].eid_; 
    size_type a_idx = edges_[i].a_.index(); //Something wrong with the edge constructor
    size_type b_idx = edges_[i].b_.index();
    //std::cout << id << ", " << a_idx << ", " << b_idx << std::endl;
    return Edge(id, a_idx, b_idx);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const { // HW0: YOUR CODE HERE
    // Make the tuples for comparison
    key_pair_type a_to_b = std::make_tuple(a.index(), b.index());
    key_pair_type b_to_a = std::make_tuple(b.index(), a.index());
    //If this a.idx and b.idx exist
    if (nodes_to_eid.find(a_to_b) != nodes_to_eid.end() || nodes_to_eid.find(b_to_a) != nodes_to_eid.end()) {
      return true;
    } else {
      return false;
    }
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
  Edge add_edge(const Node& a, const Node& b) {// HW0: YOUR CODE HERE
    key_pair_type a_to_b = std::make_tuple(a.index(), b.index());
    key_pair_type b_to_a = std::make_tuple(b.index(), a.index()); 
    if (nodes_to_eid.find(a_to_b) != nodes_to_eid.end()) {
      return Edge(nodes_to_eid[a_to_b], a.index(), b.index()); 
    } else if (nodes_to_eid.find(b_to_a) != nodes_to_eid.end()) {
      return Edge(nodes_to_eid[b_to_a], a.index(), b.index()); 
    } else { //Make a new edge
      edges_.push_back(internal_edge(eid_to_assign, a, b)); //Add a new node instance to the vector
      eid_to_idx[eid_to_assign] = num_edges()-1; //This might be a little funky
      auto temp_tup = std::make_tuple(a.index(), b.index()); 
      nodes_to_eid[temp_tup] = eid_to_assign; //Creating the map
      ++eid_to_assign; 
      return Edge(eid_to_assign - 1, a.index(), b.index());        // Invalid Edge
    }
  }
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() { // HW0: YOUR CODE HERE
    //Clear all of the edges
    edges_.clear(); //Clear the edge vector
    nodePairs_.clear();
    eid_to_idx.clear();
    //node_to_eid.clear();
    eid_to_assign = 0;

    //Clear all of the nodes
    nodes_.clear();
    assert(num_nodes() == 0); 
    assert(num_edges() == 0); 
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

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  class internal_node {
    Point point; //Figure out how to get the position
    size_type uid; // The unique identification for this node
    /** Constructor */
    internal_node(Point position, size_type uid_) {
      point = position;
      uid = uid_;
    }
    friend class Graph;
  };

 private:// HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // Internal type for nodes
  //class internal_node {
  //  Point point; //Figure out how to get the position
  //  size_type uid; // The unique identification for this node
  //  /** Constructor */
  //  internal_node(Point position, size_type uid) {
  //    point = position;
  //    uid = uid;
  //  }
  //  friend class Node;
  //};

  std::vector<internal_node> nodes_; //Create a vector of nodes
  //size_type size_;
  //size_type next_uid_;

class internal_edge {
  size_type eid_;
  Node a_; 
  Node b_;
  /** Constructor */
  internal_edge(size_type eid, Node a, Node b)
    : eid_(eid), a_(a), b_(b) {}
  friend class Graph;
};


  /** Edge attributes */
  size_type idx_;
  size_type eid_to_assign; //Edge identification number
  std::vector<internal_edge> edges_; // Vector of edges
  std::map<size_type, size_type> eid_to_idx; // eid:idx
  std::map<key_pair_type, size_type> nodes_to_eid;
//  std::unordered_map<NodeSet, size_type, HashEdge> node_to_eid; // (na, nb):eid
  std::vector<std::tuple<size_type, size_type>> nodePairs_; //Creating a vector of node pairs. Idx will correspond with edge.
 // std::unordered_map<std::tuple<Node, Node>, size_type> node_to_eid;
 // std::unordered_map<Edge, size_type> node_to_idx; // (na, nb):eid

  // Disable copy and assignment of a graph
  // I think this is the same as deleting it?
  //Graph(const Graph&) = delete;
  //Graph& operator = (const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
