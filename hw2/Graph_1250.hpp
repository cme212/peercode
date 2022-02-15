#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/Color.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
// template <typename V>
template <typename V, typename E>
class Graph
{
private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

  /** Predeclaration of node internal struct. */
  struct internal_node;
  struct adjacent_node;

  /** Predeclaration of edge helper structs. */
  struct NodePair;
  struct hash_func;

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
  using edge_value_type = E; // equivalent to typedef E edge_value_type;

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
  using size_type = unsigned; // evaluates to unsigned int

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  { // initialized using defaults of STL containers
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODE PROXY
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>
  {
  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    Graph *graph_;

    // This node's unique identification number
    size_type uid_;

    /** Private Constructor */
    Node(const Graph *graph, size_type uid)
        : graph_(const_cast<Graph *>(graph)), uid_(uid)
    {
    }

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
      graph_ = nullptr; // initialize node unattached to a graph
      uid_ = -1;        // initialize uid that is easily identified as invalid
    }

    /** Return this node's unique identifier. */
    size_type getid() const
    {
      return uid_;
    }

    /** Return internals associated with node */
    internal_node fetch()
    {
      size_type idx = graph_->uid_to_idx[uid_];
      return graph_->nodes_[idx];
    }

    /** Helper method to find initial position of node */
    Point init_position()
    {
      size_type idx = graph_->uid_to_idx[uid_];
      return graph_->nodes_[idx].p0_;
    }

    /** Return this node's position (modifiable version). */
    Point &position()
    {
      // extract uid from uid_to_idx hash map
      size_type idx = graph_->uid_to_idx[uid_];

      // extract point info from node internals
      Point &node_position = graph_->nodes_[idx].point_;

      return node_position;
    }

    /** Return this node's position (const version). */
    const Point &position() const
    {
      // extract uid from uid_to_idx hash map
      size_type idx = graph_->uid_to_idx[uid_];

      // extract point info from node internals
      Point &node_position = graph_->nodes_[idx].point_;

      return node_position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      // extract uid from uid_to_idx hash map
      size_type idx = graph_->uid_to_idx[uid_];

      // check if index valid
      size_type graph_size = graph_->uid_to_idx.size();
      if (((int)idx >= 0) and (idx < graph_size))
      {
        return idx;
      }

      return size_type(-1);
    }

    /** Return this node's value. */
    node_value_type &value()
    {
      unsigned int idx = graph_->uid_to_idx[uid_];
      internal_node &node = graph_->nodes_[idx];
      return node.val_;
    }

    const node_value_type &value() const
    {
      unsigned int idx = graph_->uid_to_idx[uid_];
      internal_node &node = graph_->nodes_[idx];
      return node.val_;
    }

    /** Return degree = number of incident edges */
    size_type degree() const
    {
      size_type idx = this->index();
      return graph_->nodes_[idx].adj_nodes_.size();
    }

    /** Start of the incident iterator */
    incident_iterator edge_begin() const
    {
      // std::cout << "hey i'm incident iterator begin()" << std::endl;
      return incident_iterator(this->graph_, *this, 0);
    }

    /** End of incident iterator */
    incident_iterator edge_end() const
    {
      // std::cout << "hey i'm incident iterator end()" << std::endl;
      unsigned int idx = graph_->uid_to_idx[uid_];
      unsigned int iter_len = graph_->nodes_[idx].adj_nodes_.size();
      return incident_iterator(this->graph_, *this, iter_len);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      // check unique id equality
      bool eq_uid = uid_ == n.getid();

      // check same graph
      bool eq_graph = this->graph_ == n.graph_;

      if (eq_uid && eq_graph)
      {
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
    bool operator<(const Node &n) const
    {
      // check unique id equality
      bool lessthan_uid = uid_ < n.getid();

      // check same graph
      bool lessthan_graph = this->graph_ < n.graph_;

      if (lessthan_uid || lessthan_graph)
      {
        return true;
      }
      return false;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    // map container affords us size() method
    return uid_to_idx.size();
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
  Node add_node(const Point &position,
                const node_value_type &val = node_value_type())
  {

    // declare and initialize uid, idx
    size_type new_uid = next_uid_;
    size_type idx = uid_to_idx.size();
    uid_to_idx[new_uid] = idx;

    // instantiate node internals and append to graph
    internal_node new_node = {position, new_uid, val};
    nodes_.push_back(new_node);

    // increment next_uid_ tracker
    ++next_uid_;

    return Node(this, new_uid);
  }

  /** Remove node from the graph, returning q if successful, else 0
   * @param[in] node The node to be remove
   * @post new num_nodes() == old num_nodes() - 1
   * @post all adjacent edges to node removed
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(num_nodes()) at most, but hopefully a lot less
   */
  size_type remove_node(const Node &n)
  {
    if (has_node(n))
    {
      // find unique id, idx
      size_type uid = n.getid();
      size_type idx = n.index();

      // extract adjacent nodes
      std::vector<unsigned int> &adj_nodes = nodes_[idx].adj_nodes_;

      // remove all edges in this node's adjacency lists
      for (auto it = adj_nodes.begin(); it != adj_nodes.end();) // no increment
      {
        // extract uid, node of adjacent edge
        size_type adj_uid = *it;
        size_type adj_idx = uid_to_idx.at(adj_uid);
        Node adj_n = Node(this, adj_uid);

        // find corresponding ref in other adjacency
        std::vector<unsigned int> &adj_nodes2 = nodes_[adj_idx].adj_nodes_;
        auto it2 = adj_nodes2.begin();
        while (it2 != adj_nodes2.end())
        {
          if (uid == *it2)
          {
            *it2 = adj_nodes2.back();
            adj_nodes2.pop_back(); // void return
            break;
          }
          else
          {
            ++it2;
          }
        }

        // remove ref from this adjacency list
        *it = adj_nodes.back();
        adj_nodes.pop_back(); // void return

        // remove adjacent edge
        size_type res = remove_edge(n, adj_n);
        (void)res; // quiet compiler warning
      }

      // swap indices in nodes_
      internal_node tmp_node = nodes_.back();
      nodes_[idx] = nodes_.back();
      nodes_.pop_back();

      // update uid_to_idx for swapped internal node
      uid_to_idx[tmp_node.uid_] = idx;

      // erase node from uid_to_idx
      uid_to_idx.erase(uid);

      return 1;
    }

    return 0;
  }

  /** Wrapper function for remove_node if passed node_iterator
   * @param[in] n_it is a node_iterator pointing to node to be removed
   * @post should not increment iterator if assert statement returns true
   */
  node_iterator remove_node(node_iterator n_it)
  {
    Node n = *n_it;
    if (remove_node(n))
    {
      return (*this).node_begin();
    }

    return ++n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {

    const size_type node_id = n.getid(); // access id from node class
    if (uid_to_idx.find(node_id) != uid_to_idx.end())
    { // search hash map
      return true;
    }

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
    size_type node_id = nodes_[i].uid_; // extract uid_
    return Node(this, node_id);         // instantiate proxy node
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
  class Edge : private totally_ordered<Edge>
  {
  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // two nodes associated with edge
    Node a_, b_;

    /** Private Constructor */
    Edge(const Node &a, const Node &b)
        : a_(a), b_(b)
    {
    }

  public:
    /** Construct an invalid Edge. */
    Edge()
        : a_(), b_()
    {
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return a_;
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return b_;
    }

    /** Return length (euclidean distance) of this edge */
    double length() const
    {
      return norm_2(a_.position() - b_.position());
    }

    /** Return or assign this edge's value */
    edge_value_type &value()
    {
      Graph *g = a_.graph_;
      size_type eid = g->nodes_to_eid[NodePair(a_, b_)];
      return g->eid_to_val[eid];
    }

    const edge_value_type &value() const
    {
      Graph *g = a_.graph_;
      size_type eid = g->nodes_to_eid[NodePair(a_, b_)];
      return g->eid_to_val[eid];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      // check equality of nodes in both directions
      bool case1 = (a_ == e.node1()) and (b_ == e.node2());
      bool case2 = (b_ == e.node1()) and (a_ == e.node2());

      if (case1 or case2)
      { // note, checks graph equality via node equality
        return true;
      }

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      // check sum of node ids
      size_type sum_ab = a_.getid() + b_.getid();
      size_type sum_ab2 = e.node1().getid() + e.node2().getid();
      bool lessthan_sum = sum_ab < sum_ab2;

      // check smallest node id
      size_type min_ab = std::min(a_.getid(), b_.getid());
      size_type min_ab2 = std::min(e.node1().getid(), e.node2().getid());
      bool lessthan_min = min_ab < min_ab2;

      // check if edges on different graphs
      bool lessthan_graph = a_.graph_ < e.node1().graph_;

      if (lessthan_sum || lessthan_min || lessthan_graph)
      {
        return true;
      }
      return false;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return eid_to_idx.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return Edge(edges_[i].node1(), edges_[i].node2());
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    // check both directions in nodes_to_eid hash map
    bool case1 = nodes_to_eid.find(NodePair(a, b)) != nodes_to_eid.end();
    bool case2 = nodes_to_eid.find(NodePair(b, a)) != nodes_to_eid.end();

    if (case1 or case2)
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
    // instantiate new edge
    Edge e(a, b);

    // check if preexisting
    if (!has_edge(a, b))
    {

      // append new edge to edges_
      edges_.push_back(e);

      // declare and initialise eid, idx
      size_type new_eid = next_eid_;
      size_type idx = eid_to_idx.size();

      eid_to_idx[new_eid] = idx;
      nodes_to_eid[NodePair(a, b)] = new_eid;
      eid_to_val[new_eid] = edge_value_type();

      // instantiate adjacent nodes (for inc_it_)
      unsigned int adj_to_a = b.uid_;
      unsigned int adj_to_b = a.uid_;
      // std::cout << "node 1: " << a.uid_ << std::endl;
      // std::cout << "node 2: " << b.uid_ << std::endl;

      // update node internals to reflect adjacents
      size_type idx_a = uid_to_idx[a.uid_];
      size_type idx_b = uid_to_idx[b.uid_];

      nodes_[idx_a].adj_nodes_.push_back(adj_to_a);
      nodes_[idx_b].adj_nodes_.push_back(adj_to_b);

      // increment next_eid, adj_nodes
      ++next_eid_;

      return e;
    }

    return e;
  }

  /** Removes an edge from our graph
   * @param[in] a of Node type, one node of edge to be removed
   * @param[in] b of Node type, one node of edge to be removed
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an unsigned int that can be evaluated as a boolean (1 for success)
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node &a, const Node &b)
  {
    if (has_edge(a, b))
    {
      // identify eid, idx
      size_type eid = nodes_to_eid[NodePair(a, b)];
      size_type idx = eid_to_idx.at(eid);

      // extract uids from node args
      size_type uid_a = a.getid();
      size_type uid_b = b.getid();

      // update adjacency list for first node
      size_type idx_a = uid_to_idx.at(uid_a);
      std::vector<unsigned int> &adj_a = nodes_[idx_a].adj_nodes_;

      for (int i = 0; i < (int)adj_a.size(); ++i)
      {
        if (adj_a[i] == uid_b)
        {
          adj_a[i] = adj_a.back();
          adj_a.pop_back(); // void return
          break;
        }
      }

      // update adjacency list for second node
      size_type idx_b = uid_to_idx.at(uid_b);
      std::vector<unsigned int> &adj_b = nodes_[idx_b].adj_nodes_;

      for (int j = 0; j < (int)adj_b.size(); ++j)
      {
        if (adj_b[j] == uid_a)
        {
          adj_b[j] = adj_b.back();
          adj_b.pop_back(); // void return
          break;
        }
      }

      // swap indices in edges_
      Edge tmp_edge = edges_.back();
      edges_[idx] = edges_.back();
      edges_.pop_back();

      // update eid_to_idx for swapped edge
      NodePair tmp_node_pair(tmp_edge.a_, tmp_edge.b_);
      eid_to_idx[nodes_to_eid[tmp_node_pair]] = idx;

      // erase from eid_to_idx map
      eid_to_idx.erase(eid);

      // erase from eid_to_val map
      eid_to_val.erase(eid);

      // erase from nodes_to_eid map
      nodes_to_eid.erase(NodePair(a, b));

      return 1;
    }

    return 0;
  }

  /** Wrapper function for remove_edge if passed edge_type
   * @param[in] e edge to be removed
   */
  size_type remove_edge(const Edge &e)
  {
    Node a = e.node1();
    Node b = e.node2();

    return remove_edge(a, b);
  }

  /** Wrapper function for remove_edge if passed
   * @param[in] e_it edge iterator pointing to edge to be removed
   * @post should not increment iterator if this function returns true
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    Edge e = *e_it;
    if (remove_edge(e))
    {
      return (*this).edge_begin();
    }
    return ++e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    nodes_.clear();
    edges_.clear();
    uid_to_idx.clear();
    eid_to_idx.clear();
    nodes_to_eid.clear();
    eid_to_val.clear();
    next_uid_ = 0;
    next_eid_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>
  {
  private:
    friend class Graph;

    Graph *graph_;
    unsigned int node_idx_;
    unsigned int nodes_len_;

    NodeIterator(const Graph *pgraph,
                 unsigned int node_idx)
        : graph_(const_cast<Graph *>(pgraph)), node_idx_(node_idx),
          nodes_len_(graph_->uid_to_idx.size())
    {
    }

  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    // Define operator++: increments to the next node in graph
    NodeIterator &operator++()
    {
      if (node_idx_ < nodes_len_)
      {
        ++node_idx_;
      }
      return *this; // return lvalue reference
    }

    // Define operator==: equality and inequality of node iterators
    bool operator==(const NodeIterator &node_iter) const
    {
      bool eq_graph_ptr = (this->graph_ == node_iter.graph_);
      bool eq_node_idx = (this->node_idx_ == node_iter.node_idx_);
      return eq_graph_ptr && eq_node_idx;
    }

    // Define operator*: dereferencing operator for node iterator
    node_type operator*() const
    {
      unsigned int node_uid = graph_->nodes_[node_idx_].uid_;
      return node_type(graph_, node_uid);
    }
  };

  // Return iterator pointing to first node of graph
  NodeIterator node_begin() const
  {
    return NodeIterator(this, 0);
  }

  // Return iterator pointing one-past-last node of graph
  NodeIterator node_end() const
  {
    return NodeIterator(this, uid_to_idx.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>
  {
  private:
    friend class Graph;
    Graph *graph_; // use graph ptr instead of const_iterator obj
    node_type root_node_;
    unsigned int iter_id_;
    unsigned int root_idx_;
    unsigned int iter_len_;

    IncidentIterator(const Graph *pgraph, node_type root_node, unsigned int iter_id)
        : graph_(const_cast<Graph *>(pgraph)), root_node_(root_node), iter_id_(iter_id),
          root_idx_(graph_->uid_to_idx[root_node_.uid_]),
          iter_len_(graph_->nodes_[root_idx_].adj_nodes_.size())
    {
    }

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
    }

    /** Define operator++: increments to the next incident node */
    IncidentIterator operator++()
    {
      // std::cout << " ++ incident iterator" << std::endl;
      if (iter_id_ < iter_len_)
      { // move to next adj node
        ++iter_id_;
      }
      return *this; // returns lvalue reference
    }

    /** Define operator++: equality of incident iterators */
    bool operator==(const IncidentIterator &incident_iter) const
    {
      // std::cout << " == incident iterator" << std::endl;
      bool eq_graph_ptr = (this->graph_ == incident_iter.graph_);
      bool eq_root_node = (this->root_node_ == incident_iter.root_node_);
      bool eq_iter_id = (this->iter_id_ == incident_iter.iter_id_);
      return eq_graph_ptr && eq_root_node && eq_iter_id;
    }

    /** Define operator*: dereferencing operator for incident iterator
     * @pre callers should use * operator with at or after end()
     * */
    edge_type operator*() const
    {
      // std::cout << " * incident iterator" << std::endl;
      unsigned int adj_node_uid = graph_->nodes_[root_idx_].adj_nodes_[iter_id_];
      return Edge(root_node_, Node(graph_, adj_node_uid));
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>
  {
  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    typename std::vector<Edge>::const_iterator edge_it_;

    EdgeIterator(typename std::vector<Edge>::const_iterator edge_it)
        : edge_it_(edge_it) {}

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

    /** Define operator*: dereferencing operator for edge iterator */
    Edge operator*() const
    {
      return Edge((*edge_it_).a_, (*edge_it_).b_);
    }

    /** Define operator++: increments to the next edge */
    EdgeIterator &operator++()
    {
      edge_it_++;
      return *this; // return lvalue reference
    }

    /** Define operator==: equality of edge iterators */
    bool operator==(const EdgeIterator &edge_iter) const
    {
      return this->edge_it_ == edge_iter.edge_it_;
    }
  };

  /** Return iterator pointing to first edge of graph */
  edge_iterator edge_begin() const
  {
    return edge_iterator(edges_.begin());
  }

  /** Return iterator pointing one-past-last edge of graph */
  edge_iterator edge_end() const
  {
    return edge_iterator(edges_.end());
  }

private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node
  {
    Point p0_;
    Point point_;
    size_type uid_;
    node_value_type val_;
    std::vector<unsigned int> adj_nodes_;

    internal_node(const Point &position,
                  size_type uid,
                  const node_value_type &val)
        : p0_(position), point_(position), uid_(uid), val_(val), adj_nodes_() {}
  };

  struct NodePair
  {
    size_type min_uid;
    size_type max_uid;

    NodePair(const Node &a, const Node &b)
    {
      min_uid = std::min(a.uid_, b.uid_);
      max_uid = std::max(a.uid_, b.uid_);
    }

    // compare keys in hash collisions
    bool operator==(const NodePair &np) const
    {
      return min_uid == np.min_uid and max_uid == np.max_uid;
    }

    // less than operator for edges
    bool operator<(const NodePair &np) const
    {
      return min_uid < np.min_uid; // compare smallest nodes
    }
  };

  struct hash_func
  {
    size_type operator()(const NodePair &np) const
    {
      // Hash function taken from boost library
      size_type seed = 0;
      std::hash<size_type> hasher;
      seed ^= hasher(np.min_uid) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= hasher(np.max_uid) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };

  struct adjacent_node
  {
    unsigned int adj_eid_; // id of edge
    unsigned int adj_uid_; // id of node
    adjacent_node(unsigned int adj_eid, unsigned int adj_uid)
        : adj_eid_(adj_eid), adj_uid_(adj_uid) {}
  };

  /** Declare containers for nodes. */
  size_type next_uid_ = 0;
  std::vector<internal_node> nodes_;
  std::unordered_map<size_type, size_type> uid_to_idx;

  /** Declare containers for edges. */
  size_type next_eid_ = 0;
  std::vector<Edge> edges_;
  std::unordered_map<size_type, size_type> eid_to_idx;
  std::unordered_map<size_type, edge_value_type> eid_to_val;
  std::unordered_map<NodePair, size_type, hash_func> nodes_to_eid;
};

#endif // CME212_GRAPH_HPP
