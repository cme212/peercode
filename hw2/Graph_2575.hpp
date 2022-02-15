#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// Add the template
template <typename V, typename E>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // Internal node representation
  struct node_rep;
  // Internal edge representation
  struct edge_rep;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  typedef V node_value_type;
  typedef E edge_value_type;

  using node_map_type = const std::unordered_map<unsigned,node_rep>*;

  using edge_map_type = const std::unordered_map<unsigned, edge_rep>*;

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
  // Store nodes in a map, edges in a map, node_edge_check for adding edge
  Graph() 
      : node_map{}, total_nodes(0), total_edges(0), node_size(0), edge_size(0), 
      edge_map{}, nid_eid_map{}, active_nodes{}, active_edges{} {
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
  class Node: private totally_ordered<Node> {
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
    // Public constructor, creates invalid node
    Node() {
    }

    /** Non-const version to return this node's position, 
     * allows for modification. */
    Point& position() {
      return graph_->node_map.at(ind_).p;
    }

    /** Return this node's position. Checks the node map at index. */
    const Point& position() const {
      return graph_->node_map.at(ind_).p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->node_map.at(ind_).user_idx;
    }

    /** Return this node's value of type V, allows for modifying value. */
    node_value_type& value() {
      return graph_->node_map.at(ind_).val;
    }

    /** Return this node's value of type V, doesn't allow for modification. */
    const node_value_type& value() const {
      return graph_->node_map.at(ind_).val;
    }

    /** Returns this node's degree, the number of adjacent nodes. */
    size_type degree() const {
      return graph_->node_map.at(ind_).adj_nodes.size();
    }

    /** Returns beginning incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this,graph_->node_map.at(ind_).
      adj_nodes.begin(),graph_->node_map.at(ind_).adj_nodes.end());
    }

    /** Returns end incident iterator 
    (corresponding to the last incident node)
    */
    incident_iterator edge_end() const {
      return IncidentIterator(this,graph_->node_map.at(ind_).adj_nodes.end(),
      graph_->node_map.at(ind_).adj_nodes.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Checks whether the graph and index are equal.
      if ((n.graph_ == graph_) && 
      (n.ind_ == ind_)) {
        return true;
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
      if ((ind_ < n.ind_) || (n.graph_ != graph_ && ind_ == n.ind_)) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the graph
    Graph* graph_;
    // Node unique is
    size_type ind_;
     /** Private Constructor */
    Node(const Graph* graph, size_type ind)
      : graph_(const_cast<Graph*>(graph)), ind_(ind) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position, the node value
   * node value is set to a default if not supplied
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
  const node_value_type& val = node_value_type()) {
    // Inserts a node to the node_map with empty adjacency set to start
    std::unordered_set<size_type> empty_adj_set;
    // Note total nodes will be unique
    node_map.insert({total_nodes,node_rep(position,total_nodes,node_size,val,
    empty_adj_set)});
    // Add to active nodes
    active_nodes.emplace_back(total_nodes);
    // Adds 1 to the node_size
    node_size += 1;
    // Add 1 to the total nodes
    total_nodes += 1;
    // Returns node
    return Node(this,total_nodes-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Searches in the map.
    if (node_map.find(n.ind_) == node_map.end()) {
      return false;
    }
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,active_nodes.at(i));
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

      return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Checks for undirected nature of the graph. Note this
      // implicity checks graph pointers
      if (((e.node1() == node1_) && (e.node2() == node2_)) ||
      ((e.node2() == node1_) && (e.node1() == node2_))) {
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
      // Tests if the two nodes are less than each other.
      if ((node1_ < e.node1()) && (node2_ < e.node2())) {
          return true;
      }
      return false;
    }

    /** Returns the edge length as the Euclidean distance between two nodes. */
    double length() const {
      return norm_2(node1_.position()-node2_.position());
    }

    /** Returns the edge value, allows for modification. */
    edge_value_type& value() {
      // Note the graphs will be equal for node 1 and node 2.
      // Get the edge id
      size_type edge_id = 
      node1_.graph_->nid_eid_map.at({node1_.ind_,node2_.ind_});
      // Get the edge value
      return node1_.graph_->edge_map.at(edge_id).edge_val;
    }

    /** Returns the egde value, does not allow for modification. */
    const edge_value_type& value() const {
      size_type edge_id = 
      node1_.graph_->nid_eid_map.at({node1_.ind_,node2_.ind_});
      return node1_.graph_->edge_map.at(edge_id).edge_val;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Define an edge as containing two nodes that point to the same Graph.
    node_type node1_;
    node_type node2_;

    /* Private constructor */
    Edge(const node_type node1, const node_type node2)
      : node1_(node1), node2_(node2) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Get the edge unique id
    size_type e_id = active_edges.at(i);
    // Get the user ids for the nodes
    size_type user_n1  = node_map.at(edge_map.at(e_id).node_ind1).user_idx;
    size_type user_n2 = node_map.at(edge_map.at(e_id).node_ind2).user_idx;
    // Return the edge
    return Edge(node(user_n1),node(user_n2));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1) because we are using unordered maps.
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Considers undirected nature of the graph
    // Check if the two nodes of this edge are added to the graph
    if ((node_map.at(a.ind_).adj_nodes.find(b.ind_) 
    == node_map.at(a.ind_).adj_nodes.end()) &&
    (node_map.at(b.ind_).adj_nodes.find(a.ind_) 
    == node_map.at(b.ind_).adj_nodes.end())) {
      return false;
    }
    return true;
  }

  /** Add an edge to the graph, return the current edge if it already exists.
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
  Edge add_edge(const Node& a, const Node& b,
  const edge_value_type& val = edge_value_type()) {
    // Check if the edge exists.
    if (has_edge(a,b)) {
      return Edge(a,b);
    }
    // Add to the edge map if it doesn't exist
    edge_map.insert({total_edges,edge_rep(a.ind_,b.ind_,val,
    total_edges,edge_size)});

    // Add to the node adjacency sets
    node_map.at(a.ind_).adj_nodes.insert(b.ind_);
    node_map.at(b.ind_).adj_nodes.insert(a.ind_);

    // Add to node id -> edge id map
    nid_eid_map.insert({{a.ind_,b.ind_},total_edges});

    // Add to the vector of active edges
    active_edges.emplace_back(total_edges);
    
    // Add one to the total edges and edges size
    total_edges += 1;
    edge_size += 1;
    return Edge(a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_map.clear();
    edge_map.clear();
    nid_eid_map.clear();
    active_nodes.clear();
    active_edges.clear();
    node_size = 0;
    edge_size = 0;
    total_edges = 0;
    total_nodes = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  //private totally_ordered<Node
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;// Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Dereferencing operator. 
    Returns the node corresponding to the graph at a given index.
    */
    Node operator*() const {
      // Note here iter_idx corresponds to the user index
      // returns a node at this index
      return Node(graph_ptr,
      graph_ptr->node_map.at(graph_ptr->active_nodes[iter_idx]).node_idx);
    }
    
    /** Advances the iterator. */
    NodeIterator& operator++() {
      if (iter_idx < graph_ptr->num_nodes()) {
        ++iter_idx;
      }
      return *this;
    }

    /** Iterators are equal if graphs and node indices are equal. */
    bool operator==(const NodeIterator& node_iter) const {
      return node_iter.graph_ptr == graph_ptr 
      && node_iter.iter_idx == iter_idx;
    }


   private:
    friend class Graph;

    const graph_type* graph_ptr;
    size_type iter_idx;

    // Private constructor
    NodeIterator(const graph_type* graph_ptr_, size_type iter_idx_) :
    graph_ptr(graph_ptr_), iter_idx{iter_idx_} {
    }
  };

  /** Begin node iterator, node at 0. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** End node iterator, one past the end. */
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes());
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
    using reference         = Edge&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Dereferencing operator, returns edge attached
    to initial node and current adjacent node.
    */
    Edge operator*() const {
      // return Edge(*node_ptr, (*node_ptr).graph_
      // ->node((*node_ptr).graph_->node_map.at(*it).user_idx));
      return Edge(*node_ptr,Node((*node_ptr).graph_,*it));
    }

    /** Increment operator. */
    IncidentIterator& operator++() {
      if (it != end) {
        ++it;
      }
      return *this;
    }

    /** Equality operator. */
    bool operator==(const IncidentIterator& inc_iter) const {
      return node_ptr == inc_iter.node_ptr 
      && it == inc_iter.it && end == inc_iter.end;
    }

   private:
    friend class Graph;

    // stores a node pointer, iterator to beginning of adjacency set
    // iterator to end of adjacency set
    const Node* node_ptr;
    std::unordered_set<size_type>::iterator it;
    std::unordered_set<size_type>::iterator end;   

    IncidentIterator(const Node* node_ptr_, 
    std::unordered_set<size_type>::iterator it_,
    std::unordered_set<size_type>::iterator end_) :
    node_ptr(node_ptr_), it{it_}, end{end_} {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to element
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy



    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Returns the two nodes connected to this edge. */
    Edge operator*() const {
      // Get the two node indices using the active edge vector
      size_type n1_ind = 
      graph_ptr->edge_map.at(graph_ptr->active_edges[iter_idx]).node_ind1;
      size_type n2_ind = 
      graph_ptr->edge_map.at(graph_ptr->active_edges[iter_idx]).node_ind2;
      // Return the edge corresponding to two nodes
      return Edge(Node(graph_ptr,n1_ind),
      Node(graph_ptr,n2_ind));
    }

    /** ++ operator. */
    EdgeIterator& operator++() {
      if (iter_idx < graph_ptr->num_edges()) {
        ++iter_idx;
      }
      return *this;
    }

    /** Equal if they have the same graph and edge index. */
    bool operator==(const EdgeIterator& edge_iter) const {
      return edge_iter.graph_ptr == graph_ptr && 
      edge_iter.iter_idx == iter_idx;
    }

   private:
    friend class Graph;

    // Edge iterator has graph pointer and iteration index.
    const graph_type* graph_ptr;
    size_type iter_idx;

    EdgeIterator(const graph_type* graph_ptr_, size_type iter_idx_) :
    graph_ptr(graph_ptr_), iter_idx{iter_idx_} {
    }
  };

  /** Edge iterator begins at first index. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /** Edge iterator end. */
  edge_iterator edge_end() const {
    return EdgeIterator(this,this->num_edges());
  }

   /** 
    * @brief Removes an edge from the graph.
    * 
    * @param[in] e An edge of the graph.
    * @return 1 if the graph has the edge and it was removed, 0 else.
    * 
    * @pre _e_ is a valid edge.
    * @pre 0 <= _e_ index < num_edges.
    * @post If old has_edge(_e_), new has_edge(_e_) == false
    * @post If old has_edge(_e_), new num_edges() == old num_edges() - 1.
    *       Else,                        new num_edges() == old num_edges().
    * @post If old has_edge(_e_), _e_ is removed from data structures
    * containing this edge.
    * 
    * Uses the swap/pop procedure to remove the specified edge from the
    * vector of active_edges, updates the corresponding data structures
    * with correct unique ids and user indexes such that
    * the updated graph does not contain the edge.
    *
    * Complexity: O(log(num_edges()) (due to searching in a map).
    */
  size_type remove_edge(const Edge& e) {
    // Checks if we have the edge, O(1)
    if (has_edge(e.node1_,e.node2_)) {
      
      // Get the edge unique id
      size_type e_ind = nid_eid_map.at({e.node1_.ind_,e.node2_.ind_});
      // Swap the user indices in the edge map
      edge_map.at(active_edges.back()).user_idx = edge_map.at(e_ind).user_idx;
      // Swap and pop to remove the specified edge from active edges
      std::swap(active_edges[edge_map.at(e_ind).user_idx],active_edges.back());
      active_edges.pop_back();
      
      // Erase from edge_map, n_id_eid map
      // Map erase by iterator is average case constant time
      edge_map.erase(edge_map.find(e_ind));
      nid_eid_map.erase(nid_eid_map.find({e.node1_.ind_,e.node2_.ind_}));
      // Note unordered set erase is average constant time in number of
      // elements removed
      node_map.at(e.node1_.ind_).adj_nodes.erase(e.node2_.ind_);
      node_map.at(e.node2_.ind_).adj_nodes.erase(e.node1_.ind_);

      // Subtract 1 from edge size
      edge_size -= 1;
      return 1;
    }
    // else return 0
    return 0;
  }

  /** 
    * @brief Removes an edge from the graph given two nodes of the edge.
    * 
    * @param[in] n1 node1 of the edge.
    * @param[in] n2 node2 of the edge.
    * @return 1 if the graph has the edge and it was removed, 0 else.
    * 
    * @pre _n1_ and _n2_ are valid nodes
    * @pre edge(_n1_,_n2_) is a valid edge
    * @pre 0 <= _n1_ index < num_nodes()
    * @pre 0 <= _n2_ index < num_nodes()
    * @post If old has_edge(_n1_,_n2_), new has_edge(_n1_,_n2_) == false
    * @post If old has_edge(_n1_,_n2_), new num_edges() == old num_edges() - 1.
    *       Else,                        new num_edges() == old num_edges().
    * @post If old has_edge(_n1_,_n2_), edge(_n1_,_n2_) is removed from 
    * data structures containing this edge.
    * 
    * Uses the swap/pop procedure to remove the specified edge from the
    * vector of active_edges, updates the corresponding data structures
    * with correct unique ids and user indexes such that
    * the updated graph does not contain the edge.
    *
    * Complexity: O(log(num_edges()) (due to searching in a map).
    */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // check if the edge exists
    if (has_edge(n1,n2)) {
      // Get the unique edge is
      size_type e_ind = nid_eid_map.at({n1.ind_,n2.ind_});
      // Swap the user indices in the edge map
      edge_map.at(active_edges.back()).user_idx = edge_map.at(e_ind).user_idx;
      // Swap and pop to remove the specified edge from active edges
      std::swap(active_edges[edge_map.at(e_ind).user_idx],active_edges.back());
      active_edges.pop_back();

      // Erase from all data structures containing this edge
      edge_map.erase(edge_map.find(e_ind));
      nid_eid_map.erase(nid_eid_map.find({n1.ind_,n2.ind_}));
      node_map.at(n2.ind_).adj_nodes.erase(n1.ind_);
      node_map.at(n1.ind_).adj_nodes.erase(n2.ind_);
      
      // Subtract 1 from edge size
      edge_size -= 1;
      return 1;
    }
    return 0;
  }

  /** 
    * @brief Removes an edge from the graph given an edge iterator.
    * 
    * @param[in,out] e_it edge iterator.
    * @return 1 if the graph has the edge and it was removed, 0 else.
    * 
    * @pre _e_it_ is a valid edge iterator
    * @post If old has_edge(_e_it_), new has_edge(_e_it_) == false
    * @post If old has_edge(*_e_it_), new num_edges() == old num_edges() - 1.
    *       Else,                        new num_edges() == old num_edges().
    * @post If old has_edge(*_e_it), edge(*_e_it_) is removed from 
    * data structures containing this edge.
    * @post *_e_it_ either points to the beginning edge or is unchanged.
    * 
    * Uses previous remove_edge methods to remove the edge pointed
    * to by an edge_iterator.
    *
    * Complexity: O(log(num_edges())
    */
  edge_iterator remove_edge(edge_iterator e_it) {
    // Remove edge
    size_type removed = remove_edge(*e_it);
    // Successful, return beginning
    if (removed != 0) {
      return this->edge_begin();
    }
    // Return iterator if unsuccessful
    // eg the graph doesn't have the edge
    return e_it;
  }

  /** 
    * @brief Removes a node from the graph given the node.
    * 
    * @param[in] n node to remove.
    * @return 1 if the graph has the edge and it was removed, 0 else.
    * 
    * @pre _n_ is a valid node.
    * @pre 0 <= _n_ index() < num_nodes()
    * @post If old has_node(_n_), new has_node(_n_) == false
    * @post If old has_node(_n_), new num_edges() == old num_edges() - 1.
    *       Else,                        new num_edges() == old num_edges().
    * @post If old has_node(_n_), node(_n_) is removed from 
    * data structures containing this edge.
    * @post Adjacent edges to _n_ are removed if old has_node(_n_)
    * 
    * Uses the swap/pop procedure to remove the specified node from the
    * vector of active_node, updates the corresponding data structures
    * with correct unique ids and user indexes such that
    * the updated graph does not contain the edge. Also uses
    * remove_edge to remove all adjacent edges. 
    *
    * Complexity: O(log(num_edges()) + _n_.degree())
    */
  size_type remove_node(const Node& n) {
    // Check if the graph contains the node.
    if (has_node(n)) {
      // Get the adjacent nodes, copy is O(degree)
      // Wasn't able to use iterator because elements
      // in iterator are modified in loop
      // Remove n.ind_ from any adjacent nodes in remove edge
      auto adj_node_copy = node_map.at(n.ind_).adj_nodes;
      // Removes adjacent edges
      for (auto adj_node: adj_node_copy) {
        remove_edge(n,node(node_map.at(adj_node).user_idx));
      }
      // Swaps user indices
      node_map.at(active_nodes.back()).user_idx = node_map.at(n.ind_).user_idx;
      // Swap and pop to remove this node from list of active nodes
      std::swap(active_nodes[node_map.at(n.ind_).user_idx]
      ,active_nodes.back());
      active_nodes.pop_back();
      // Erases from node map
      node_map.erase(node_map.find(n.ind_)); 
      // Adjust node size
      node_size -= 1;
      return 1;
    }
    return 0;  
  }

  /** 
    * @brief Removes a node from the graph given a node iterator.
    * 
    * @param[in,out] n_it node iterator.
    * @return 1 if the graph has the node and it was removed, 0 else.
    * 
    * @pre _n_it_ is a valid node iterator
    * @post If old has_node(*_n_it_), new has_node(*_n_it_) == false
    * @post If old has_node(*_n_it_), new num_nodes() == old num_nodes() - 1.
    *       Else,                        new num_nodes() == old num_nodes().
    * @post If old has_node(*_n_it), node(*_n_it_) is removed from 
    * data structures containing this edge.
    * @post *_n_it_ either points to the beginning edge or is unchanged.
    * @post Adjacent to edges to *n_it are removed if old has_node(*n_it)
    * 
    * Uses previous remove_node methods to remove the node pointed
    * to by a node_iterator.
    *
    * Complexity: O(log(num_edges())
    */
  node_iterator remove_node(node_iterator node_it) {
    // Remove node
    size_type removed = remove_node(*node_it);
    // Check if removed
    if (removed != 0) {
      return this->node_begin();
    }
    return node_it;
  }



 private:

  // structure for internal nodes, contains a point and node index
  // value and set of adjacent nodes
  struct node_rep {
    Point p;
    size_type node_idx; // Unique node index
    size_type user_idx; // Outward facing user index
    node_value_type val; // Node value
    std::unordered_set<size_type> adj_nodes; // Adjacent nodes
    node_rep(Point curr_p, size_type curr_idx, size_type uid, 
    node_value_type v, 
    std::unordered_set<size_type> adj_nodes_) : p(curr_p),
    node_idx(curr_idx), user_idx(uid), val(v), adj_nodes(adj_nodes_) {}
  };

  // structure for internal edge, contains two node indices
  struct edge_rep {
    size_type node_ind1; // Node1 index
    size_type node_ind2; // Node2 index
    edge_value_type edge_val; // Edge value
    size_type edge_idx; // Unique edge index
    size_type user_idx; // Outward facing user index
    edge_rep(size_type i1, size_type i2, edge_value_type v,
    size_type eid, size_type uid) : node_ind1(i1), node_ind2(i2)
    , edge_val(v), edge_idx(eid), user_idx(uid) {}
  };

  // Node unique index to node rep
  std::unordered_map<size_type,node_rep> node_map; 
  size_type total_nodes;
  size_type total_edges;
  size_type node_size;
  size_type edge_size;
  // Edge unique index to edge rep
  std::unordered_map<size_type, edge_rep> edge_map;
  // {n1,n2} to edge unique index
  std::map<std::set<size_type>, size_type> nid_eid_map;
  // Active nodes and edges that can be accessed
  std::vector<size_type> active_nodes;
  std::vector<size_type> active_edges;
};

#endif // CME212_GRAPH_HPP
