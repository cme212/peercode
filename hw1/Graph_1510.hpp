#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int>
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

  using node_value_type = V;
  
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
  using pair_type = std::pair<size_type, size_type>;
  using vec_pair_type = std::vector<std::pair<size_type, size_type>>;
  using set_type = std::set<size_type>;
  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), edge_map_(), incident_map_(),
      next_node_id_(0), next_edge_id_(0) {}

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
      graph_ = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      //if invalid node, no position
      assert(graph_ != nullptr);
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
      //if invalid node, no index
      assert(graph_ != nullptr);
      return node_id_; // node_id_ is set to be the index
    }

    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the value of the node. if invalid, raises assert error.
     */
    node_value_type& value() {
      //if invalid node, no value
      assert(graph_!=nullptr);
      
      for (size_type i=0; i < graph_->size(); i++) {
        if (graph_->nodes_[i].node_id == node_id_) {
          return graph_->nodes_[i].val;
        }
      }
      assert(false);
    }
    
    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the value of the node as a const. if invalid, raises assert error.
     */
    const node_value_type& value() const {
      //if invalid node, no value
      assert(graph_ != nullptr);
      
      for (size_type i=0; i < graph_->size(); i++) {
        if (graph_->nodes_[i].node_id == node_id_) {
          return graph_->nodes_[i].val;
        }
      }
      assert(false);
    };
    
    /** Return the value of a node.
     * @pre current node is a valid node in this graph.
     * @return the number of edges for this node size.
     * @post 0 <= return <= num_nodes() - 1
     */
    size_type degree() const {
      //if invalid node or valid node doesn't have any edges
      if (!graph_ or
          (graph_->incident_map_.find(node_id_)==graph_->incident_map_.end())){
        return 0;
      }
      //if node is in incident_map
      auto& v = graph_->incident_map_.at(node_id_);
      return v.size();
    }
    
    /** Return an incident iterator pointing to the first edge of the node.
     * @pre current node is a valid node in this graph.
     * @return an incident iterator pointing to this graph's this node and its first edge.
     * @return return is an invalid incident iterator if a node doesn't have any edges.
     */
    IncidentIterator edge_begin() const {
      //if invalid node, then no edges nor incident iterator
      assert(graph_ != nullptr);
      
      if (graph_->incident_map_.find(node_id_) == graph_->incident_map_.end()) {
        return IncidentIterator();
      }
      //if node is in incident_map
      auto& v = graph_->incident_map_.at(node_id_);
      return IncidentIterator(graph_, node_id_, v.begin());
    }
    
    /** Return an incident iterator pointing past the last edge of the node.
     * @pre current node is a valid node in this graph.
     * @return an incident iterator pointing to this graph's this node
     *         and past its last edge.
     * @return return is an invalid incident iterator if a node doesn't have
     *         any edges.
     */
    IncidentIterator edge_end() const {
      //if invalid node, then no edges nor incident iterator
      assert(graph_ != nullptr);
      
      //node does not have any edges
      if (graph_->incident_map_.find(node_id_) == graph_->incident_map_.end()) {
        return IncidentIterator();
      }
      //if node is in incident_map
      auto& v = graph_->incident_map_.at(node_id_);
      return IncidentIterator(graph_, node_id_, v.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (!graph_ and !n.graph_) {
        return true;
      }
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
      if (!graph_ and !n.graph_) {
        return false;
      }
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
  Node add_node (const Point& position,
                 const node_value_type& val = node_value_type()) {
    internal_node intNode {position, next_node_id_, val};
    nodes_.push_back(intNode);
    //add to incident map with empty vector
    if (incident_map_.find(next_node_id_) == incident_map_.end()) {
      vec_pair_type vec {};
      std::pair<size_type, vec_pair_type> kvPair {next_node_id_, vec};
      incident_map_.insert(kvPair);
    }
    
    ++next_node_id_;
    return Node(this, next_node_id_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_; //check if n points to this Graph
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, first_node_id_);
    }
    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, second_node_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (!graph_ and !e.graph_) {
        return true;
      }
      // same two nodes
      return graph_ == e.graph_ &&
              ((node1()==e.node1() && node2()==e.node2()) ||
              (node1()==e.node2() && node2()==e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (!graph_ and !e.graph_) {
        return false;
      }
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
    //This edge's nodes specified in order.
    size_type first_node_id_;
    size_type second_node_id_;
    //Private Constructor
    Edge(const Graph* graph,
         size_type edge_id,
         size_type first_node_id,
         size_type second_node_id)
         : graph_(const_cast<Graph*> (graph)),
           edge_id_(edge_id),
           first_node_id_(first_node_id),
           second_node_id_(second_node_id) {}
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
    if (i < this->edges_.size()) {
      internal_edge intEdge = this->edges_[i];
      return Edge(this, i, intEdge.nodeA_id, intEdge.nodeB_id);
    }
    return Edge(); //return invalid edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    set_type set_ {a.node_id_, b.node_id_};
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
    set_type set_ {a.node_id_, b.node_id_};
    //CASE A: Edge{a,b} already exists
    if (has_edge(a, b)) {
      size_type eid = edge_map_.at(set_);
      return Edge(this, eid, a.index(), b.index());
    }
    //CASE B: new edge
    //add edge to edges_ vector
    internal_edge intEdge {a.node_id_, b.node_id_, next_edge_id_};
    edges_.push_back(intEdge);
    //add edge to edge_map_
    std::pair<set_type, size_type> edgePair (set_, next_edge_id_);
    edge_map_.insert(edgePair);
    
    //add nodes to incident_map_
    //Case 1: Node A is in the map already; append to vec
    if (incident_map_.find(a.node_id_) != incident_map_.end()) {
      auto& a_vec = incident_map_.at(a.node_id_);
      pair_type newPair {b.node_id_, next_edge_id_};
      a_vec.push_back(newPair);
    } else { //Case 2: create Node A in map
      //pair to put into vector; create 1-element vector
      pair_type newPair{b.node_id_, next_edge_id_};
      vec_pair_type vec {newPair};
      //create pair to insert into map
      std::pair<size_type, vec_pair_type> kvPair {a.node_id_, vec};
      incident_map_.insert(kvPair);
    }
    //Case 3: Node B is in the map already
    if (incident_map_.find(b.node_id_) != incident_map_.end()) {
      auto& b_vec = incident_map_.at(b.node_id_);
      pair_type newPair2 {a.node_id_, next_edge_id_};
      b_vec.push_back(newPair2);
    } else { //Case 4: create Node B in map
      pair_type newPair {a.node_id_, next_edge_id_};
      vec_pair_type vec {newPair};
      std::pair<size_type, vec_pair_type> kvPair {b.node_id_, vec};
      incident_map_.insert(kvPair);
    }
    
    ++next_edge_id_;
    return Edge(this, next_edge_id_-1, a.index(), b.index());
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
    incident_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered <NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy
    typedef typename std::vector<internal_node>::const_iterator
                                                 vector_iterator_ ;
    
    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_ = nullptr;
    }

    /** Return a node pointed to by this NodeIterator.
     * @pre node iterator points to a valid internal node in nodes_ .
     * @return a Node with the node_id specified by the internal node.
     */
    Node operator*() const {
      //if this iterator is invalid, raise error.
      assert(graph_ != nullptr);
      internal_node intNode = *nit_;
      return graph_->node(intNode.node_id);
    }
    /** Return this NodeIterator by reference.
     * @pre node iterator points to a valid internal node in nodes_
     *      in this graph, not past it.
     * @return this NodeIterator by reference, incremented from pre-call.
     */
    NodeIterator& operator++() {
      //if this iterator is invalid, raise error.
      assert(graph_ != nullptr);
      nit_++;
      return *this;
    }
    /** Return bool indicating if this NodeIterator equals other.
     * @pre NodeIterators are valid.
     * @return 1 if they point to the same graph and element in nodes_,
     *         0 otherwise.
     *         1 if they're both invalid iterators.
     */
    bool operator==(const NodeIterator& other) const {
      if (!graph_ and !other.graph_) {
        return true;
      }
      return nit_ == other.nit_ and graph_ == other.graph_;
    }

   private:
    friend class Graph;
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //Iterator over nodes_ vector for graph_
    vector_iterator_ nit_;
    NodeIterator(const Graph* graph, vector_iterator_ nit)
    : graph_(const_cast<Graph*> (graph)), nit_(nit) {}
  };

  /** Return a NodeIterator pointing to the first elem of nodes_.
   * @pre nodes_ is a valid vector. its size must be >= 0.
   * @return a NodeIterator pointing to the beginning of nodes_.
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }
  /** Return a NodeIterator pointing past the last elem of nodes_.
   * @pre nodes_ is a valid vector. its size must be >= 0.
   * @return a NodeIterator pointing to the end of nodes_.
   */
  NodeIterator node_end() const {
    return NodeIterator(this, nodes_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    
    typedef typename std::vector<std::pair<size_type, size_type>>
                        ::const_iterator vector_iterator_ ;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graph_ = nullptr;
    }

    /** Return an Edge pointed to by this IncidentIterator.
     * @pre edge iterator points to a valid internal edge in this graph.
     * @pre egde iterator is valid.
     * @return an Edge with the edge ID from the internal edge it pts to.
     * @post returned edge must return this nodeID_ as its 1st node.
     */
    Edge operator*() const {
      //confirm this iterator is valid.
      assert(graph_ != nullptr);
      
      if (pairs_ == graph_->incident_map_.at(nodeID_).end()) {
        return Edge();
      }
      size_type second_nodeID_ = (*pairs_).first;
      size_type eid = (*pairs_).second;
      return Edge(graph_, eid, nodeID_, second_nodeID_);
    }
    
    /** Return this after incrementing.
     * @pre edge iterator is valid.
     * @return return == this with pairs_ incremented.
     */
    IncidentIterator& operator++() {
      //confirm this iterator is valid.
      assert(graph_ != nullptr);
      
      pairs_++;
      return *this;
    }
    
    /** Return a boolean indicating if this and other iterators are equal.
     * @return 1 if they point to the same graph and element in edges_,
     *         0 otherwise.
     *         1 if they're both invalid iterators.
     */
    bool operator==(const IncidentIterator& other) const {
      if (!graph_ and !other.graph_) {
        return true;
      }
      return graph_ == other.graph_
             && nodeID_ == other.nodeID_
             && pairs_ == other.pairs_;
    }

   private:
    friend class Graph;
    //Pointer to the Graph this Edge belongs to
    Graph* graph_;
    //nodeID whose edges this iterator goes over.
    int nodeID_;
    //iterator for the edges.
    vector_iterator_ pairs_;
    
    IncidentIterator(const Graph* graph, int nID, vector_iterator_ p)
    : graph_(const_cast<Graph*> (graph)), nodeID_(nID), pairs_(p) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered <EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    typedef typename std::vector<internal_edge>::const_iterator
                                              vector_iterator_ ;
    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_ = nullptr;
    }

    /** Return an edge this is pointing to.
     * @pre this is a valid iterator.
     * @return an edge pointed to by this iterator.
     * @post return==a valid edge.
     */
    Edge operator*() const {
      //confirm that this is a valid iterator.
      assert(graph_ != nullptr);
      internal_edge intEdge = *eit_;
      return Edge(graph_, intEdge.edge_id,
                  intEdge.nodeA_id, intEdge.nodeB_id);
    }
    /** Return this.
     * @pre this is a valid iterator.
     * @return this after incrementing @a eit_.
     * @post return is a valid iterator.
     */
    EdgeIterator& operator++() {
      //confirm that this is a valid iterator.
      assert(graph_ != nullptr);
      eit_++;
      return *this;
    }
    /** Return a boolean indicating if 2 iterators are equal..
     * @return 1 if 2 iterators' graph_ and eit_  are equal.
     * @return 0 if they're not.
     * @return if both this and other are invalid, return 1.
     */
    bool operator==(const EdgeIterator& other) const {
      //if both iterators are invalid, return true
      if (!graph_ and !other.graph_) {
        return true;
      }
      return graph_ == other.graph_
             and eit_ == other.eit_;
    }

   private:
    friend class Graph;
    //Pointer to the Graph this Node belongs to
    Graph* graph_;
    //Iterator pointing to an int. edge in edges_
    vector_iterator_ eit_;
    
    EdgeIterator(const Graph* graph, vector_iterator_ eit)
    : graph_(const_cast<Graph*> (graph)), eit_(eit) {}
  };

  /** Return an EdgeIterator pointing at the 1st elem of edges_.
   * @pre edges_ is a valid vector. its size must be >= 0.
   * @return an EdgeIterator pointing to the 1st elem of edges_.
   * @post a valid EdgeIterator
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }
  /** Return an EdgeIterator pointing at the last elem of edges_.
   * @pre edges_ is a valid vector. its size must be >= 0.
   * @return an EdgeIterator pointing to the last elem of edges_.
   * @post a valid EdgeIterator
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:
  struct internal_node {
    Point point; //the Point the node stores
    size_type node_id; //unique identification of the node
    node_value_type val;
  };
  struct internal_edge {
    size_type nodeA_id;
    size_type nodeB_id;
    size_type edge_id; //unique identification of the edge
  };
  
  //stores all nodes as internal nodes
  std::vector<internal_node> nodes_;
  //stores all edges as internal edges
  std::vector<internal_edge> edges_;
  //stores key=node pairs, val=edge ID
  std::map<set_type, size_type> edge_map_;
  //stores key=node, val=vec of (node2, edgeID) pairs
  std::unordered_map<size_type, vec_pair_type> incident_map_;

  size_type next_node_id_;
  size_type next_edge_id_;
};

#endif // CME212_GRAPH_HPP
