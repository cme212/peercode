#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
 
/** @file Graph.hpp
 * @brief An undirected graph type
 */   

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int> // default is int, if V is not specified
class Graph {
public:
  using node_value_type = V;

private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node; // Stores point, index in vector, and uid

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

  /** Construct an empty graph. Its graph UID will be 0 */
  Graph()
    // HW0: YOUR CODE HERE
    // empty maps (both the ordered and the unordered), since initially there are no nodes and no edges
    // empty map uid->idx, since initially there are no nodes
    // next node UID is 0 (for the very first node)
    // no nodes yet, so nodes_ is empty.
    : nodes_(), next_node_uid_(0), node_uid_to_idx_(), edge_map_unordered_(), edge_map_ordered_() {}

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
    Node() 
      // invalid node points to no graph. Its UID is irrelevant so just set it to 0
      : uid_(0), graph_(nullptr) {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[graph_->node_uid_to_idx_.at(uid_)].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->node_uid_to_idx_.at(uid_);
    }

    

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const; 

    /**
     * METHOD: value()
     * (non-const version)
     * 
     * @pre this node is a valid node in some graph
     * @return a reference to the internal node's value 
     *         that you can use to actually change that value
     */
    node_value_type& value() {
      return const_cast<graph_type*>(graph_)->nodes_[graph_->node_uid_to_idx_.at(uid_)].val_;
    }

    /** METHOD: value()
     * (const version)
     * 
     * @pre this node is a valid node in some graph
     * @return a reference to the interla node's value
     *         for read-purposes only. You cannot use it to change the internal value.
     */
    const node_value_type& value() const {
      return graph_->nodes_[graph_->node_uid_to_idx_.at(uid_)].val_;
    }

    /**
     * METHOD: degree()
     * 
     * @return 0 if node is invalid, else return degree in the graph
     */
    size_type degree() const {
      if (graph_==nullptr) {return 0;}
      return graph_->edge_map_unordered_.at(uid_).size();
    }

    /**
     * METHOD: edge_begin()
     * 
     * @pre node is a valid node in some graph
     * 
     * @return an incident iterator to start iterating over that node's adjacent edges
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    /**
     * METHOD: edge_end()
     * 
     * @pre node is a valid node in some graph
     * @return the incident iterator pointing to the end of that node's adjacent edges.
     *         It doesn't point to any edge in particular. (Points one past the end).
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_,
                               uid_,
                               graph_->edge_map_unordered_.at(uid_).size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_!=n.graph_) {
        return false;
      }
      if (graph_ == nullptr) {
        if (n.graph_ == nullptr) {return uid_ == n.uid_;}
        return false;
      }
      if (n.graph_ == nullptr) {return false;}
      if (graph_->node_uid_to_idx_.find(uid_) == graph_->node_uid_to_idx_.end() ||
          n.graph_->node_uid_to_idx_.find(n.uid_) == n.graph_->node_uid_to_idx_.end()) {
        return false;
      }
      return graph_->node_uid_to_idx_.at(uid_)==n.graph_->node_uid_to_idx_.at(n.uid_);
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
      if (graph_==n.graph_) {
        return uid_<n.uid_;
      }
      return graph_<n.graph_;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    friend class Graph::Edge;
    friend class Graph::NodeIterator;
    friend class Graph::IncidentIterator;

    // Construct a node pointing to a given graph with a given uid
    Node(const graph_type* graph_pointer, size_type node_uid) : uid_(node_uid) , graph_(graph_pointer) {}

    size_type uid_;
    const graph_type* graph_; // Node is allowed to modify its parent graph, through value()
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return (size_type) nodes_.size();
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
    // HW0: YOUR CODE HERE
    // Create a new node with a given position.
    // Optional value. If not given, value defaults to the default of whatever type it is (e.g. 0 for int).
    Node new_node(this, next_node_uid_);
    internal_node new_internal_node(next_node_uid_, position, size(), val); // pass val by value, not reference
    node_uid_to_idx_[next_node_uid_] = size(); 
    nodes_.push_back(new_internal_node);
    // That vector currently has no neighbours, so its list of neighbours is an empty vector.
    edge_map_unordered_[next_node_uid_] = std::vector<size_type>();
    edge_map_ordered_[next_node_uid_] = std::vector<size_type>();
    next_node_uid_++;
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_==this) &&
           (node_uid_to_idx_.find(n.uid_) != node_uid_to_idx_.end());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    Node new_node(this, nodes_[i].uid_);
    return new_node;
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
    Edge() 
      // HW0: YOUR CODE HERE
      : node1_(0), node2_(0), graph_(nullptr) {}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {return false;}
      if (this->node1_ == e.node1_ && this->node2_ == e.node2_) {return true;}
      if (this->node1_ == e.node2_ && this->node2_ == e.node1_) {return true;}
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // If in same graph, compare the endpoint nodes' UIDs lexicographically.
      if (this->graph_==e.graph_) {
        if (this->node1_ < this->node2_) {
          if (e.node1_ < e.node2_) {
            if (this->node1_ < e.node1_) {return true;}
            else if (this->node1_ == e.node1_) {return this->node2_ < e.node2_;}
            else {return false;}
          }
          else {
            if (this->node1_ < e.node2_) {return true;}
            else if (this->node1_ == e.node2_) {return this->node2_ < e.node1_;}
            else {return false;}
          }
        }
        else {
          if (e.node1_ < e.node2_) {
            if (this->node2_ < e.node1_) {return true;}
            else if (this->node2_ == e.node1_) {return this->node1_ < e.node2_;}
            else {return false;}
          }
          else {
            if (this->node2_ < e.node2_) {return true;}
            else if (this->node2_ == e.node2_) {return this->node1_ < e.node1_;}
            else {return false;}
          }
        }
      }
      // If not in same graph, compare the graphs.
      return this->graph_<e.graph_;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    friend class Graph::IncidentIterator;
    friend class Graph::EdgeIterator;

    Edge(const graph_type* graph_pointer, size_type node1_uid, size_type node2_uid)
      : node1_(node1_uid), node2_(node2_uid), graph_(graph_pointer) {}

    size_type node1_;
    size_type node2_;
    const graph_type* graph_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type count = 0;
    for(auto iter = edge_map_ordered_.begin();
        iter != edge_map_ordered_.end();
        ++iter) 
    {
      count += (iter->second).size();
    }
    return count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type count = 0;
    for(auto iter = edge_map_ordered_.begin();
        iter != edge_map_ordered_.end();
        ++iter) 
    {
      if (count + (iter->second).size() > i) {
        return Edge(this, iter->first, (iter->second)[i-count]);
      }
      count += (iter->second).size();
    }
    //Only reach this line if i > num_edges()
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type i=0; i<edge_map_unordered_.at(a.uid_).size(); i++) {
      if (edge_map_unordered_.at(a.uid_)[i] == b.uid_) return true;
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
    Edge new_edge(this, a.uid_, b.uid_);
    for (size_type i=0; i<edge_map_unordered_.at(a.uid_).size(); i++) {
      if (edge_map_unordered_.at(a.uid_)[i] == b.uid_) return new_edge;
    }

    //Edge not present in graph
    edge_map_unordered_.at(a.uid_).push_back(b.uid_);
    edge_map_unordered_.at(b.uid_).push_back(a.uid_);
    if (a.uid_ <= b.uid_) edge_map_ordered_.at(a.uid_).push_back(b.uid_);
    if (b.uid_ <= a.uid_) edge_map_ordered_.at(b.uid_).push_back(a.uid_);

    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    node_uid_to_idx_.clear();
    next_node_uid_=0;
    edge_map_unordered_.clear();
    edge_map_ordered_.clear();
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
    NodeIterator() : index_(0), graph_(nullptr) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * METHOD: operator*()
     * 
     * @pre nodeiterator points to existing graph
     * @pre nodeiterator is not pointing at the end (i.e. one past the last node),
     *      but rather it is pointing to a valid node in the graph.
     *      So, for instance, *(graph.node_end()) is illegal.
     * @post the method does not modify *this
     * 
     * @return The node that the iterator is pointing to.
     */
    Node operator*() const {
      return graph_->node(index_);
    }

    /**
     * METHOD: operator++()
     * 
     * @pre nodeiterator points to existing graph
     * 
     * @post nodeiterator now points to the next node,
     *       or to the end() if it was one before the end,
     *       or it stays at the end() if it was already at the end().
     * 
     * @return a reference to this same nodeiterator, after the ++ has been performed.
     */
    NodeIterator& operator++() {
      if (index_ < graph_->num_nodes()) {index_++;}
      NodeIterator& ni = *this;
      return ni;
    }

    /**
     * METHOD: operator==()
     * 
     * @param[in] @a ni is a reference to a NodeIterator that the function
     *            does not modify
     * @return whether or not the two nodeiterators point to the same node 
     *         (or are both at the end) of the same graph.
     * @post the method does not modify *this
     */
    bool operator==(const NodeIterator& ni) const {
      return (graph_ == ni.graph_) && (index_ == ni.index_);
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // Nodeiterator is just a pointer to its parent graph
    // and the index of the node it's currently pointing to.
    // When it's at the end, index_ == graph_->size()
    // and the nodeiterator doesn't point to any node in particular
    size_type index_;
    const graph_type* graph_;

    // Constructor for node iterator given index of node it's pointing to, and the graph it's in.
    NodeIterator(size_type index, const graph_type* graph_pointer) : index_(index), graph_(graph_pointer) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * METHOD: node_begin()
   * 
   * @return a nodeiterator of the graph pointing to the first node.
   *         or to the end() in the case that the graph has zero nodes
   * 
   * @post the method does not modify *this
   */
  node_iterator node_begin() const {
    return node_iterator(0, this);
  }

  /**
   * METHOD: node_end()
   * 
   * @return a nodeiterator of the graph pointing one past the last node,
   *         (so *(graph.node_end()) is illegal)
   * 
   * @post the method does not modify *this
   */
  node_iterator node_end() const {
    return node_iterator(num_nodes(), this);
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    //There's no source node and no index node, so the values of these attribtues don't matter
    //The graph pointer must point nowhere.
    IncidentIterator() : graph_(nullptr), source_(0), index_(0) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * METHOD: operator*()
     * 
     * @pre incidentiterator belongs to an existing graph
     *                       and points to a valid edge (i.e. is not at the end())
     * 
     * @return the Edge it points to
     * 
     * @post this method does not modify @this
     */
    Edge operator*() const {
      if (graph_==nullptr) return Edge();
      if (index_ < graph_->edge_map_unordered_.at(source_).size()) {
        return Edge(graph_, source_, graph_->edge_map_unordered_.at(source_)[index_]);
      }
      return Edge();
    }

    /**
     * METHOD: operator++()
     * 
     * @post the iterator stays where it's at, if it was at the end(),
     *       or if it was an invalid iterator. Otherwise, it now points
     *       to the next available edge, or to the end() if it was previously
     *       at the last available edge.
     * 
     * @return a reference to *this
     */
    IncidentIterator& operator++() {
      if (graph_==nullptr) {
        IncidentIterator& ii2 = *this;
        return ii2;
      }
      if (index_ < graph_->edge_map_unordered_.at(source_).size()) index_++;
      IncidentIterator& ii2 = *this;
      return ii2;
    }

    /**
     * METHOD: operator==()
     * 
     * @param[in] @a ii2 another incidentiterator to compare this to.
     * 
     * @post this method does not modidfy this
     * 
     * @return whether or not the iterators point to the same edge
     *         of the same source node in the same graph. Note that it is not
     *         enough for the edges to match; the "root node" (source_) must match as well.
     *         Two invalid incidentiterators will also evaluate to being equal.
     */
    bool operator==(const IncidentIterator& ii2) const {
      return (graph_ == ii2.graph_) && 
             (source_ == ii2.source_) &&
             (index_ == ii2.index_);
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    friend class Graph::Node;

    // Construct incidentiterator given the graph, the source node, and the index of that
    // source node's list of adjacent edges
    IncidentIterator(const graph_type* graph_pointer, size_type source, size_type index) :
                     graph_(graph_pointer), source_(source), index_(index) {}

    const graph_type* graph_;
    size_type source_; // source node's UID
    size_type index_; // index over the vector of neighbours of the source node
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    //Points to no graph. the map_iterator_ is a default map iterator. The index defaults to 0
    EdgeIterator() : graph_(nullptr), map_iterator_(), index_(0) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * METHOD: operator*()
     * 
     * @pre The edgeiterator belongs to a valid graph
     * 
     * @post this method does not modify *this
     * 
     * @return the Edge this edgeiterator points to in its parent graph, or an invalid Edge
     *         object in the case that the edgeiterator is pointing at the end.
     *         E.G. if g is an empty graph, then *(g.edge_begin()) is an invalid Edge object.
     *         E.G. For any graph g, *(g.edge_end()) is an invalid edge object.
     */
    Edge operator*() const {
      // If we're at end, return invalid edge
      if (map_iterator_ == graph_->edge_map_ordered_.end()) return Edge();
      // If index is after the end of the source node's vector of neighbours, return invalid edge
      // (However, this case should never happen) in practice
      if (index_ >= map_iterator_->second.size()) return Edge();

      // Get node endpoints
      return Edge(graph_, map_iterator_->first, map_iterator_->second[index_]);
    }

    /**
     * METHOD: operator++()
     * 
     * @post The edgeiterator is not modified in the case that it is an invalid edgeiterator object
     *       or in the case that it was already pointing at the end.
     *       Otherwise, the edgeiterator now points to the next available edge, or to the end in the
     *       case that it was previously pointing to the last available edge.
     * 
     * @return a reference to *this
     */
    EdgeIterator& operator++() {
      // If invalid edgeiterator object, return
      if (graph_==nullptr) {
        EdgeIterator& ei2 = *this;
        return ei2;
      }
      // If we're at end, return 
      if (map_iterator_ == graph_->edge_map_ordered_.end()) {
        EdgeIterator& ei2 = *this;
        return ei2;
      }

      // If index is not at end of source node's vector of neighbours, just increment index
      if (index_ < map_iterator_->second.size()-1) {
        index_++;
        EdgeIterator& ei2 = *this;
        return ei2;
      }

      // Find next node with nonzero amount of neighbours 
      ++map_iterator_;
      while (map_iterator_ != graph_->edge_map_ordered_.end()) {
        if (map_iterator_->second.size()==0) {
          ++map_iterator_;
          continue;
        }
        index_ = 0;
        EdgeIterator& ei2 = *this;
        return ei2;
      }

      // If we're here, then we found no more edges. Return the end()
      map_iterator_ = graph_->edge_map_ordered_.end();
      index_ = 0;
      EdgeIterator& ei2 = *this;
      return ei2;
    }

    /**
     * METHOD: operator==()
     * 
     * @param[in] @a ei1, a edgeiterator object to compare this to.
     * 
     * @post this method does not modify *this
     * 
     * @return whether or not the edgeiterators point to the same edge in the same graph.
     *         Behaviour when one or both of the edgeiterator objects is invalid is undefined.
     */
    bool operator==(const EdgeIterator& ei2) const {
      return (graph_ == ei2.graph_) &&
             (map_iterator_ == ei2.map_iterator_) &&
             (index_ == ei2.index_);
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    // Initialize an edge iterator pointing to the start, given its parent graph.
    // The graph_pointer parameter must be a valid graph.
    EdgeIterator(const graph_type* graph_pointer) : graph_(graph_pointer),
                                                    map_iterator_(graph_pointer->edge_map_ordered_.begin()),
                                                    index_(0) {}

    const graph_type* graph_;
    std::map<size_type, std::vector<size_type>>::const_iterator map_iterator_; // over edge_map_ordered_
    size_type index_; // iterate over index in the vector values of the map
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * METHOD: edge_begin()
   * @pre this is a valid graph object. (It can be empty and/or edgeless).
   * 
   * @post the method does not modify *this
   * 
   * @return an edge iterator pointing to the beginning of the graph's edges.
   *         So *(graph.edge_begin()) must always be a valid edge, EXCEPT WHEN
   *         the graph is edgeless, in which case edge_begin()==edge_end().
   */
  edge_iterator edge_begin() const {
    EdgeIterator ei(this);

    // Find first node with nonzero amount of neighbours
    while (ei.map_iterator_ != edge_map_ordered_.end()) {
      if (ei.map_iterator_->second.size() > 0) {
        // Found the first edge!
        ei.index_ = 0;
        return ei;
      }
      // Keep looking...
      ++(ei.map_iterator_);
    }

    // If we reach here, then graph has no edges and ei.map_iterator_ is at edge_map_ordered_.end()
    ei.index_ = 0;
    return ei;
  }

  /**
   * METHOD: edge_end()
   * @pre This is a valid graph object.
   * 
   * @post the method does not modify *this
   * 
   * @return an edge iterator point to the end ("one past the last") of the graph's edges.
   *         So *(graph.edge_end()) is always an invalid Edge object.
   */
  edge_iterator edge_end() const {
    EdgeIterator ei(this);
    ei.map_iterator_ = edge_map_ordered_.end();
    ei.index_ = 0;
    return ei;
  }

private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    const size_type uid_;
    size_type idx_;
    const Point position_;
    node_value_type val_;

    internal_node(size_type uid, const Point& position, size_type idx, node_value_type val)
      : uid_(uid), idx_(idx), position_(position), val_(val) {}

  };

  std::vector<internal_node> nodes_;
  size_type next_node_uid_;
  std::unordered_map<size_type, size_type> node_uid_to_idx_;
  std::map<size_type, std::vector<size_type>> edge_map_unordered_; // includes (a,b) and (b,a)
  std::map<size_type, std::vector<size_type>> edge_map_ordered_; // only includes (a,b) if a<=b (the nodes' UIDs)
};

#endif // CME212_GRAPH_HPP
