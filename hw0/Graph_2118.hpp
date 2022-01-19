#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:
  //declare internal types for Graph
  struct internal_node;
  struct internal_edge;

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

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), next_node_id_(0), next_edge_id_(0)
  {
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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      //find the corresponding internal node in nodes_
      for (size_type i=0; i < graph_->size(); i++) {
        if (graph_->nodes_[i].node_id == node_id_) {
          return graph_->nodes_[i].point;
        }
      }
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // node_id_ is set to be the index
      return node_id_;
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
      // same graph, same index
      return (graph_ == n.graph_ &&
              node_id_ == n.node_id_);
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
      // Comparison based on node_id_
      if (node_id_ != n.node_id_) {
        return node_id_ < n.node_id_;
      }
      else {
        return graph_ < n.graph_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    //Pointer to the Graph this Node belongs to
    Graph* graph_;
    //This Node's unique identifier
    size_type node_id_;
    //Private Constructor
    Node(const Graph* graph, size_type node_id)
    : graph_(const_cast<Graph*> (graph)), node_id_(node_id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // size of nodes_ in this graph
    return nodes_.size();
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
  Node add_node(const Point& position) {
    // check post condition: save initial size
    internal_node intNode;
    intNode.point = position;
    intNode.node_id = next_node_id_;
    nodes_.push_back(intNode);

    ++next_node_id_;
    return Node(this, next_node_id_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //check if n points to this Graph
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Find internal edge w/ same ID in edges_, return A
      for (size_type i = 0; i < graph_->edges_.size(); i++) {
        if (edge_id_ == graph_->edges_[i].edge_id) {
          return Node(graph_, graph_->edges_[i].nodeA_id);
        }
      }
      return Node(); // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Find internal edge w/ same ID in edges_, return B
      for (size_type i = 0; i < graph_->edges_.size(); i++) {
        if (edge_id_ == graph_->edges_[i].edge_id) {
          return Node(graph_, graph_->edges_[i].nodeB_id);
        }
      }
      return Node(); // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // same two nodes
      return ((node1()==e.node1() && node2()==e.node2()) ||
              (node1()==e.node2() && node2()==e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Compared based on edge_id_
      if (edge_id_ != e.edge_id_) {
        return edge_id_ < e.edge_id_;
      }
      else {
        return graph_ < e.graph_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //This edge's unique identifier
    size_type edge_id_;
    //Private Constructor
    Edge(const Graph* graph, size_type edge_id)
    : graph_(const_cast<Graph*> (graph)), edge_id_(edge_id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // size of edges_ in this graph
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Search if edge_map_ contains a and b
    std::set<size_type> set_ {a.node_id_, b.node_id_};
    return edge_map_.find(set_) != edge_map_.end();
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
    //create a set of a and b's node ID's
    std::set<size_type> set_ {a.node_id_, b.node_id_};
    //the edge a-b already exists
    if (has_edge(a, b)) {
      return Edge(this, edge_map_.at(set_));
    }
    
    //else: add edge to graph
    internal_edge intEdge;
    intEdge.nodeA_id = a.node_id_;
    intEdge.nodeB_id = b.node_id_;
    intEdge.edge_id = next_edge_id_;
    edges_.push_back(intEdge);

    std::pair<std::set<size_type>, size_type> edgePair (set_, next_edge_id_);
    edge_map_.insert(edgePair);
    
    ++next_edge_id_;
    return Edge(this, next_edge_id_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset all private members of this Graph
    next_node_id_ = 0;
    next_edge_id_ = 0;
    nodes_.clear();
    edges_.clear();
    edge_map_.clear();
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

 private:

  struct internal_node {
    Point point; //the Point the node stores
    size_type node_id; //unique identification of the node
  };
  
  struct internal_edge {
    size_type nodeA_id;
    size_type nodeB_id;
    size_type edge_id; //unique identification of the edge
  };
  
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::map<std::set<size_type>, size_type> edge_map_;
  size_type next_node_id_;
  size_type next_edge_id_;
};

#endif // CME212_GRAPH_HPP
