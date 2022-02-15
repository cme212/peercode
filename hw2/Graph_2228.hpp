#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Predeclare the internal struct
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the value in node */
  typedef V node_value_type;
  /** Type of the value in edge */
  typedef E edge_value_type;

  /** Type of the value in node */
  //using node_value_type = V;

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
  Graph() : internodes_(), interedges_(), edge_nodes_() {}

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
      // Creating an invalid node
      graphptr_ = nullptr;
      id_ = 0;
    }

    /** Return this node's position. */
    Point& position() {
      return graphptr_->internodes_.at(id_).pointptr_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graphptr_->internodes_.at(id_).pointptr_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      auto it = std::find(graphptr_->node_indices_.begin(),\
       graphptr_->node_indices_.end(), id_);
      return std::distance(graphptr_->node_indices_.begin(), it);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Returns the value attribute of this node */
    node_value_type& value() {
      return graphptr_->internodes_.at(id_).value_;
    }

    /** Returns the value attribute of this node */
    const node_value_type& value() const {
      return graphptr_->internodes_.at(id_).value_;
    }

    // return the number of incident edges.
    size_type degree() const {
      // assert that the node isn't an invalid Node (graphptr_ == nullptr)
      assert(graphptr_);
      return graphptr_->internodes_.at(id_).incident_edges_.size();
    }

    // Start of the incident iterator.
    incident_iterator edge_begin() const {
      return IncidentIterator(graphptr_, 0, id_);
    }

    // End of incident iterator.
    incident_iterator edge_end() const {
      return IncidentIterator(graphptr_, degree(), id_);
    }

    Graph* graphptr() {
      return graphptr_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graphptr_ == n.graphptr_) && (index() == n.index());
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
      return (index() < n.index()); 
    }



   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // pointer back to the point object of this node

    // The pointer to the graph object this node belongs to
    Graph* graphptr_;
    // The ID associated with this node
    size_type id_;

    /* Private Constructor */
    Node(const Graph* graph, size_type id) : 
          graphptr_(const_cast<Graph*>(graph)), id_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return internodes_.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // The ID for the new node is set to the current number of nodes 
    size_type id = size();
    // Create a new internal node object for this node
    internodes_.insert(std::pair<size_type,internal_node>(id, 
                       internal_node(position, val)));
    node_indices_.push_back(id);
    return Node(this, id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1). 
   */
  bool has_node(const Node& n) const {
    return (n.graphptr_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // return an invalid node if i == size()
    if (i == size()) {
      return Node();
    }
    // Make sure i is a valid index
    assert(i < size());
    // Find the ith internal node within the unordered map internodes_
    return Node(this, node_indices_.at(i));
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
      // Creating an invalid edge
      graphptr_ = nullptr;
      edge_id_ = 0;
      node_a_id_ = 0;
      node_b_id_ = 0;
    }

    /** Returns the Euclidean distance between both nodes in this edge */ 
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Returns the value of this edge */ 
    edge_value_type& value() {
      return graphptr_->interedges_.at(edge_id_).value_;
    }

    /** Returns the value of this edge */ 
    const edge_value_type& value() const {
      return graphptr_->interedges_.at(edge_id_).value_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Find the node ID
      return Node(graphptr_, node_a_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Find the node ID
      return Node(graphptr_, node_b_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Make sure both edges belong to the same graph
      if (graphptr_ != e.graphptr_) { return false; }
      return ((node1() == e.node1() && node2() == e.node2()) or
              (node1() == e.node2() && node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Compare edge IDs
      if (graphptr_ != e.graphptr_) { return graphptr_ < e.graphptr_; }
      return edge_id_ < e.edge_id_;
    }

    /** Return this edge's index, a number in the range [0, graph_size). */
    size_type index() const {
      auto it = std::find(graphptr_->edge_indices_.begin(),\
       graphptr_->edge_indices_.end(), edge_id_);
      return std::distance(graphptr_->edge_indices_.begin(), it);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer to the graph object this edge belongs to 
    Graph* graphptr_;
    // ID of this edge
    size_type edge_id_;
    // ID of node a of this edge
    size_type node_a_id_;
    // ID of node b of this edge
    size_type node_b_id_;


    // Constructor
    Edge(const Graph* graph, size_type id_e, size_type id_a, size_type id_b) : 
         graphptr_(const_cast<Graph*>(graph)), edge_id_(id_e), 
         node_a_id_(id_a), node_b_id_(id_b) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return interedges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Make sure i is a valid index
    assert(i < num_edges());
    // Find the internal edge at index i
    size_type edge_id = edge_indices_.at(i);
    auto inny_edge = interedges_.at(edge_id);
    return Edge(this, edge_id, inny_edge.id_a_, inny_edge.id_b_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Seach edge_nodes_ to find whether an edge exists for nodes a and b
    std::set<size_type> node_ids = {a.id_, b.id_};
    return edge_nodes_.count(node_ids);
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
    std::set<size_type> node_ids = {a.id_, b.id_};
    // If the edge already exists
    if (has_edge(a,b)) {
      // Return the current edge
      size_type edge_id = edge_nodes_.at(node_ids);
      internal_edge* inny_edge = &interedges_.at(edge_id);
      // Switch the nodes a and b
      if (inny_edge->id_a_ == b.id_ and inny_edge->id_b_ == a.id_) { 
        inny_edge->id_a_ = a.id_;
        inny_edge->id_b_ = b.id_;
      }
      return Edge(this, edge_id, a.id_, b.id_);
    }
    // Create a new edge
    else {
      // The ID for the new edge is equal to the current size of interedges_
      size_type newid = interedges_.size();
      // Add incident edge to nodes a and b
      internodes_.at(a.id_).add_inc_edge(newid);
      internodes_.at(b.id_).add_inc_edge(newid);

      internal_edge newedge = internal_edge(a.id_, b.id_);
      // Add the new internal edge to interedges_
      interedges_.insert(std::pair<size_type, internal_edge>(newid, newedge));
      // Add the node pair and edge ID to edge_nodes_
      edge_nodes_.insert({node_ids, newid});
      edge_indices_.push_back(newid);
      return Edge(this, newid, a.id_, b.id_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internodes_.clear();
    interedges_.clear();
    edge_nodes_.clear();
    node_indices_.clear();
    edge_indices_.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graphptr_ = nullptr;
      nodeidx_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference the iterator and return a Node */
    Node operator*() const {
      return graphptr_->node(nodeidx_);
    }

    /** Increment and return the iterator */
    NodeIterator& operator++() {
      nodeidx_ += 1;
      return *this;
    }

    /** Check for iterator equality */
    bool operator==(const NodeIterator& node_iter) const {
      return (graphptr_ == node_iter.graphptr_ and 
              nodeidx_ == node_iter.nodeidx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // Graph pointer object back to the associated graph
    Graph* graphptr_;
    // node index ranging from 0 to num_nodes()
    size_type nodeidx_;

    /** Private constructor for NodeIterator */
    NodeIterator(const Graph* graphptr, size_type nodeidx) : 
      graphptr_{const_cast<Graph*>(graphptr)}, nodeidx_(nodeidx) {}

  }; 

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Returns a NodeIterator starting with the Node at index 0. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Returns a NodeIterator with index==num_nodes() (one past the end). */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
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
      graphptr_ = nullptr;
      edgeidx_ = 0;
      nodeid_ = 0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference the iterator and return an Edge */
    Edge operator*() const {
      size_type edge_id = graphptr_->internodes_.at(nodeid_).incident_edges_.at(edgeidx_);
      // if this node is node2() of this edge, return an Edge with the nodes flipped
      if (graphptr_->interedges_.at(edge_id).id_b_ == nodeid_) {
        return Edge(graphptr_, edge_id, nodeid_, graphptr_->interedges_.at(edge_id).id_a_);
      }
      // return an Edge with the original node order
      else {
        return Edge(graphptr_, edge_id, graphptr_->interedges_.at(edge_id).id_a_, 
                    graphptr_->interedges_.at(edge_id).id_b_);
      }
    }
    
    /** Increment and return the iterator */
    IncidentIterator& operator++() {
      edgeidx_ += 1;
      return *this;
    }

    /** Check for iterator equality */
    bool operator==(const IncidentIterator& inc_it) const {
      return (graphptr_ == inc_it.graphptr_ and edgeidx_ == inc_it.edgeidx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Graph pointer to associated graph object of this iterator
    Graph* graphptr_;
    // the index within incident_edges for this node. Range [0, degree())
    size_type edgeidx_;
    // The ID of the corresponding node
    size_type nodeid_;

    /** Private constructor for IncidentIterator */
    IncidentIterator(const Graph* graphptr, size_type edgeidx, size_type nodeid) : 
      graphptr_{const_cast<Graph*>(graphptr)}, edgeidx_(edgeidx), nodeid_(nodeid) {}
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
      graphptr_ = nullptr;
      edgeidx_ = 0;
    }

    /** Dereference the iterator and return an Edge */
    Edge operator*() const {
      return graphptr_->edge(edgeidx_);
    }

    /** Increment and return the iterator */
    EdgeIterator& operator++() {
      edgeidx_ += 1;
      return *this;
    }

    /** Check for iterator equality */
    bool operator==(const EdgeIterator& edge_iter) const {
      return (graphptr_ == edge_iter.graphptr_ and 
              edgeidx_ == edge_iter.edgeidx_);
    }

    private:
      friend class Graph;
      // HW1 #2: YOUR CODE HERE

      // Graph pointer to associated graph object of this iterator
      Graph* graphptr_;
      // the index of the edge within range [0, num_edges())
      size_type edgeidx_;

      /** Private constructor for EdgeIterator */
      EdgeIterator(const Graph* graphptr, size_type edgeidx) : 
        graphptr_{const_cast<Graph*>(graphptr)}, edgeidx_(edgeidx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Start of the EdgeIterator
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  // End of the EdgeIterator
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /** Remove a node from the graph, if the node exists. 
   * @return 0 if the node does not exist, 1 if the removal was successful.
   * @post graph.node(@a n.index()) == false
   *
   * First, removes all incident edges to this node using remove_edge(). 
   * Then, erases correspoding values from internodes_ and node_indices.
   *
   * Complexity: O(n.degree()), which is no more than O(num_nodes()). 
   * Removing all incident edges is O(n.degree()). Assuming the graph 
   * is sparse, n.degree() is less than num_nodes(). internodes_.erase() 
   * is O(1) since it is an unordered map. Removing from node_indices is 
   * also O(1) using the swap and pop method. 
   */
  size_type remove_node(const Node& n) {
    //check if node exists
    if (internodes_.count(n.id_) == 0) { return 0; }
    //remove all edges attached to this node
    while (n.degree() != 0) {
      auto it = n.edge_begin();
      remove_edge(*it);
    }
    //remove node
    internodes_.erase(n.id_);
    // remove node from node_indices_ using swap and pop
    auto it = node_indices_.begin() + n.index();
    *it = node_indices_.back();
    node_indices_.pop_back();
    return 1;
  }


  /** Remove a node from the graph, if the node exists. 
   * @return graph.node_begin(), a valid Node Iterator
   * @post graph.node(@a (*n_it).index()) == false
   *
   * Calls remove_node(*n_it). 
   *
   * Complexity: O((*n_it).degree()), which is no more than O(num_nodes()). 
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
  }

  /** Remove an edge from the graph, if the edge exists. 
   * @return 0 if the edge does not exist, 1 if the removal was successful.
   * @post graph.edge(e.index()) == false
   *
   * Removes the edge from all corresponding data structures: edge_nodes, 
   * interedges_, and edge_indices. 
   *
   * Complexity: O(max(a.degree(), b.degree())), which is no more than 
   * O(num_nodes() + num_edges()). Removing from edge_nodes_ and interedges_ 
   * are both O(1) since they are unordered maps. Removing from the vector 
   * incident_edges_ in remove_inc_edge is O(n.degree()) to first find the 
   * corresponding element, and then O(1) to remove it using swap and pop. 
   * Finally, removing from edge_indices_ is O(1) using swap and pop. 
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // check if this edge exists
    if (!has_edge(a,b)) { return 0; }

    // find the edge id
    std::set<size_type> node_ids = {a.id_, b.id_};
    size_type edge_id = edge_nodes_.at(node_ids);

    // erase edge from unordered_maps
    edge_nodes_.erase(node_ids);
    interedges_.erase(edge_id);

    // erase edge from node a and node b's incident edges
    internodes_.at(a.id_).remove_inc_edge(edge_id);
    internodes_.at(b.id_).remove_inc_edge(edge_id);

    // erase edge from edge_indices_ using swap and pop
    size_type edge_idx = Edge(this, edge_id, a.id_, b.id_).index();
    auto it = edge_indices_.begin() + edge_idx;
    *it = edge_indices_.back();
    edge_indices_.pop_back();
    return 1;
  }

  /** Remove an edge from the graph, if the edge exists. 
   * @return 0 if the edge does not exist, 1 if the removal was successful.
   * @post graph.edge(e.index()) == false
   *
   * Removes the edge from the graph calling remove_edge(e.node1(), e.node2())
   *
   * Complexity: O(max(a.degree(), b.degree())), which is no more than 
   * O(num_nodes() + num_edges()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph, if the edge exists. 
   * @return graph.edge_begin(), a valid Edge Iterator
   * @post graph.edge(*e_it.index()) == false
   *
   * Removes the edge from the graph calling remove_edge(*e_it)
   *
   * Complexity: O(max(a.degree(), b.degree())), which is no more than 
   * O(num_nodes() + num_edges()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** @class Graph::internal_node
   * @brief Class representing the heavyweight objects of the Node class.
   *
   * The internal_node object contains a pointer to the Point object the Node
   * belongs to. 
   */
  struct internal_node {
    // The position of this node
    Point pointptr_;
    // Vector containing IDs of all edges connected to this node
    std::vector<size_type> incident_edges_;
    // Value of this node
    node_value_type value_;
    
    // Constructor
    internal_node(const Point point, node_value_type val) : pointptr_(point), value_(val) {}

    // add an edge to incident_edges
    void add_inc_edge(size_type edge_id) {
      incident_edges_.push_back(edge_id);

    }

    void remove_inc_edge(size_type edge_id) {
      auto it = std::find(incident_edges_.begin(), incident_edges_.end(), edge_id);
      *it = incident_edges_.back();
      incident_edges_.pop_back();
    }
  };

  /** @class Graph::internal_edge
   * @brief Class representing the heavyweight objects of the Edge class.
   *
   * The internal_edge object contains a the ID values for each of the Nodes
   * the edge is attached to.
   */
  struct internal_edge {
    // The ID's of node1() and node2() of this edge
    size_type id_a_;
    size_type id_b_;

    // Value of this edge
    edge_value_type value_;

    // Constructor
    internal_edge(size_type nodea_id, size_type nodeb_id) : 
                  id_a_(nodea_id), id_b_(nodeb_id) {}
  };

  // Unordered map containing the node ID as keys and internal nodes as values
  std::unordered_map<size_type, internal_node> internodes_;
  
  // Unordered map containing the edge ID as keys and internal edges as values
  std::unordered_map<size_type, internal_edge> interedges_;
  // Map containing a set of node IDs as keys and edge IDs as values
  std::map<std::set<size_type>, size_type> edge_nodes_;

  // Vector containing the indices of each node
  std::vector<size_type> node_indices_;

  // Vector containing the indices of each edge
  std::vector<size_type> edge_indices_;
};

#endif // CME212_GRAPH_HPP