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
#include <utility> // for pair and make_pair

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
public:
  typedef V node_value_type;
  typedef E edge_value_type;

private:

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

  /**
   * Annoying technicality.
   * The IncidentIterator class uses a pointer to the internal
   * edge_map_ord_values_, which is a map from size_type to a
   * vector of pairs (size_type, edge_value_type).
   * The type of this iterator is templated because it depends
   * on edge_value_type. Thus, we define it here as a typedef.
   * msvpsei = map (from) size-type (to) vector (of) pairs
   *           (size-type, edge-value-type) iterator.
   */
  typedef typename std::map<size_type,
                            std::vector<
                            std::pair<
                            size_type,
                            edge_value_type>>>::const_iterator msvpsei_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** 
   * CONSTRUCTOR: Graph()
   * 
   * Construct an empty graph.
   * Has no nodes, the next node's UID will be 0,
   * and has no edges.
   */
  Graph() : nodes_(),
            next_node_uid_(0),
            node_uid_to_idx_(),
            edge_map_unord_(),
            edge_map_ord_values_() {}

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
    /** 
     * CONSTRUCTOR: Node
     * Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later
     */
    Node() : graph_(nullptr), uid_(0) {}


    /**
     * METHOD: position() [non-const version]
     * 
     * @return a reference to this node's position.
     */
    Point& position() {
      return const_cast<graph_type*>(graph_)->nodes_[
             graph_->node_uid_to_idx_.at(uid_)].position_;
    }

    /** 
     * METHOD: position() [const version]
     * 
     * @return a reference this node's position
     *         for read-only purposes, cannot modify it.
     */
    const Point& position() const {
      return graph_->nodes_[graph_->node_uid_to_idx_.at(uid_)].position_;
    }

    /** 
     * METHOD: index()
     * 
     * @a return this node's index,
     *    a number in the range [0, graph_size)
     */
    size_type index() const {
      return graph_->node_uid_to_idx_.at(uid_);
    }

    /**
     * METHOD: value() [non-const version]
     * 
     * @pre this node is a valid node in some graph
     * @return a reference to the internal node's value 
     *         that you can use to actually change that value
     */
    node_value_type& value() {
      return const_cast<graph_type*>(graph_)->nodes_[
             graph_->node_uid_to_idx_.at(uid_)].val_;
    }

    /** 
     * METHOD: value() [const version]
     * 
     * @pre this node is a valid node in some graph
     * @return a reference to the internal node's value
     *         for read-purposes only. 
     *         You cannot use it to change the internal value.
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
      if (graph_==nullptr) { return 0; }
      return graph_->edge_map_unord_.at(uid_).size();
    }

    /**
     * METHOD: edge_begin()
     * 
     * @pre node is a valid node in some graph
     * 
     * @return an incident iterator to start iterating over
     *         that node's adjacent edges
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    /**
     * METHOD: edge_end()
     * 
     * @pre node is a valid node in some graph
     * @return the incident iterator pointing to the end of
     *         that node's adjacent edges.
     *         It doesn't point to any edge in particular.
     *         (Points one past the end).
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_,
                               uid_,
                               graph_->edge_map_unord_.at(uid_).size());
    }

    /** 
     * METHOD: operator==()
     * Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Deal with invalid nodes first
      if (graph_!=n.graph_) { return false; }
      if (graph_ == nullptr) {
        if (n.graph_ == nullptr) { return uid_ == n.uid_; }
        return false;
      }
      if (n.graph_ == nullptr) {return false;}
      if (graph_->node_uid_to_idx_.find(uid_) ==
          graph_->node_uid_to_idx_.end() ||
          n.graph_->node_uid_to_idx_.find(n.uid_) ==
          n.graph_->node_uid_to_idx_.end()) {
        return false;
      }
      // Now deal with valid nodes
      return graph_->node_uid_to_idx_.at(uid_) ==
             n.graph_->node_uid_to_idx_.at(n.uid_);
    }

    /** 
     * METHOD: operator less than
     * 
     * Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers.
     * It need not have any geometric meaning.
     */  
    bool operator<(const Node& n) const {
      if (graph_==n.graph_) { return uid_<n.uid_; }
      return graph_<n.graph_;
    } 

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class Graph::Edge;
    friend class Graph::NodeIterator;
    friend class Graph::IncidentIterator; 

    /**
     * CONSTRUCTOR: Node()
     * 
     * Construct a node given its parent graph and UID
     */
    Node(const graph_type* graph_pointer,size_type node_uid)
            : graph_(graph_pointer) , uid_(node_uid) {}

    // Pointer to parent graph
    const graph_type* graph_;

    // Node's UID
    size_type uid_;
  };

  /** 
   * METHOD: size()
   * 
   * @return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type) nodes_.size();
  }

  /** 
   * METHOD: num_nodes()
   * Synonym for size()
   */
  size_type num_nodes() const {
    return size();
  }



  /** 
   * METHOD: add_node()
   * 
   * Add a node to the graph, returning the added node.
   * 
   * @param[in] position The new node's position
   * @paream[in] val The new node's value.
   *                 Defaults to the default value of node_value_type
   * 
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post the new node is not == any previously existing nodes.
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& val = node_value_type()) { 
    // Create new proxy node
    Node new_node(this, next_node_uid_);
    // Create new interla node.
    // Pass position and val by copy, not reference
    internal_node new_internal_node(next_node_uid_, position, size(), val); 

    node_uid_to_idx_[next_node_uid_] = num_nodes(); 
    nodes_.push_back(new_internal_node);

    // That node currently has no neighbours, so 
    // its list of neighbours is an empty vector.
    edge_map_unord_[next_node_uid_] = std::vector<size_type>();
    edge_map_ord_values_[next_node_uid_] = 
            std::vector<std::pair<size_type, edge_value_type>>();

    // Increment next node's UID
    next_node_uid_++;
    return new_node;
  }



  /** 
   * METHOD: remove_node() [given the node]
   * 
   * Remove a node from the graph.
   * 
   * @param[in] n, the node.
   * @pre n.graph_ points to this graph
   * @post !has_node(n)
   * @post If old has_node(n), then
   *       new num_nodes() == old num_nodes() -1
   *       else,
   *       new num_nodes() == old num_nodes()
   * 
   * INVALIDATION SPECIFICATIONS:
   * @post Invalidates any incidentiterator with this node
   *       as the source node.
   * @post Doesn't invalidate incidentiterators with this node
   *       as the target node, but it will change which edge
   *       they point to -- potentially changes them to
   *       source_node.edge_end()
   * @post [!!!] Might invalidate nodeiterators even if they weren't
   *       pointing to this specific node. [!!!]
   * @post Invalidates any edgeiterators that pointed to an edge with
   *       this node as an endpoint.
   * @post Doesn't invalidate edgeiterators whose edge didn't have
   *       this node as an endpoint, but it might change which edge 
   *       they point to.
   * @post If old n == old n2, then new !has_node(n2).
   * @post Invalidates any edges that had n as an endpoint.
   * @post Doesn't invalidate any other Node objects.
   * 
   * @return 1 if old has_node(n), else 0
   * 
   * Complexity: O(degree(n)+sum_i degree(n_i)),
   *             where we sum over all neighbours n_i to i
   */
  size_type remove_node(const Node& n) {
    // Check if graph has node
    if (!has_node(n)) { return 0; }

    /* Erase edges.
       Removing them as we go doesn't work because
       the IncidentIterator relies on the edges being static
       as we iterate through all neighbours. */
    std::vector<Edge> edges_to_remove;
    for (auto it = n.edge_begin();
              it != n.edge_end();
              ++it) {
      edges_to_remove.push_back(*it);
    }
    for (Edge e : edges_to_remove) { remove_edge(e); }

    // Erase from internal nodes via swap and pop.
    size_type orig_index = n.index();
    size_type orig_last_node_uid = nodes_[num_nodes()-1].uid_;
    nodes_[orig_index] = nodes_[num_nodes()-1];
    nodes_[orig_index].idx_ = orig_index;
    nodes_.pop_back();

    // Erase from node UID to idx map
    node_uid_to_idx_.at(orig_last_node_uid) = orig_index;
    node_uid_to_idx_.erase(n.uid_); // erase by key

    // Erase from ordered and unordered edge maps
    edge_map_unord_.at(n.uid_).clear();
    edge_map_unord_.erase(n.uid_);
    edge_map_ord_values_.at(n.uid_).clear();
    edge_map_ord_values_.erase(n.uid_);

    // Return 1 since 1 node was removed
    return 1;
  }

  /**
   * METHOD: remove_node() [given an iterator to it]
   * 
   * Remove a node from the graph.
   * 
   * @param[in] n_it, the node iterator.
   * @pre (*n_it).graph_ points to this graph
   * @post !has_node(*n_it)
   * 
   * @post  If old has_node(*n_it), then
   *        new num_nodes() == old num_nodes()-1
   *        else,
   *        new num_nodes() == old num_nodes()
   * 
   * INVALIDATION SPECIFICATIONS:
   * Same as remove_node(const Node&);
   * 
   * @return 1 if old has_node(*n_it), else 0
   * 
   * Complexity: O(degree(n)+sum_i degree(n_i)),
   *             where we sum over all neighbours n_i to i
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
  }





  /** 
   * METHOD: has_node()
   * 
   * Determine if a Node belongs to this Graph
   * 
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_==this) &&
           (node_uid_to_idx_.find(n.uid_) != node_uid_to_idx_.end());
  }

  /** 
   * METHOD: node()
   * 
   * @return the node with index @a i.
   * 
   * @pre 0 <= @a i < num_nodes()
   * 
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, nodes_[i].uid_);
  }





  //
  // EDGES
  // 

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. 
   * Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
  public:
    /**
     * CONSTRUCTOR: Edge
     * 
     * Construct an invalid edge. Points to no edge in no graph
     */
    Edge() : graph_(nullptr), node1_(0), node2_(0) {}

    /** 
     * METHOD: node1()
     * 
     * @return The node1 of this edge
     */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** 
     * METHOD: node2()
     * 
     * @return The node2 of this edge
     */
    Node node2() const {
      return Node(graph_, node2_);
    }

    /** 
     * METHOD: length()
     * 
     * @return the Euclidean length of the edge in 3D space
     */
    double length() const {
      return norm_2(node1().position() - node2().position());
    } 





    /**
     * METHOD: value() [non-const version]
     * 
     * @pre Valid existent edge in its graph
     * 
     * @return a reference to the edge's value (edge_value_type)
     * 
     * Complexity: O(max degree of node in graph)
     */
    edge_value_type& value() {
      // Endpoint nodes' UIDs in order. n1<=n2
      size_type n1 = (node1_ <= node2_) ? node1_ : node2_;
      size_type n2 = (node1_ <= node2_) ? node2_ : node1_;

      // Length of vector of n1's larger neighbours
      size_type l = const_cast<graph_type*>(
                    graph_)->edge_map_ord_values_.at(n1).size();

      for (size_type index = 0; index < l; index++) {
        std::pair<size_type, edge_value_type>& pair =
                const_cast<graph_type*>(
                graph_)->edge_map_ord_values_.at(n1)[index];
        if (pair.first == n2) { return pair.second; }
      }

      // If reach this line, then the edge was removed from parent graph
      // or never present. So this line is just to silence compiler warnings
      assert(false);
      return const_cast<graph_type*>(
             graph_)->edge_map_ord_values_.at(n1)[0].second;
    }

    /**
     * METHOD: value() [const version]
     * 
     * @pre Valid existent edge in its graph
     * 
     * @return a reference to the edge's value (edge_value_type)
     * 
     * Complexity: O(max degree of node in graph)
     */
    const edge_value_type& value() const {
      // Endpoint nodes' UIDs in order. n1<=n2
      size_type n1 = (node1_ <= node2_) ? node1_ : node2_;
      size_type n2 = (node1_ <= node2_) ? node2_ : node1_;

      // Length of vector of n1's larger neighbours
      size_type l = graph_->edge_map_ord_values_.at(n1).size();

      for (size_type index = 0; index < l; index++) {
        auto pair = graph_->edge_map_ord_values_.at(n1)[index];
        if (pair.first == n2) { return pair.second; }
      }

      // If reach this line, then the edge was removed from parent graph
      // or never present
      assert(false);
      return graph_->edge_map_ord_values_.at(n1)[0].second;
    }




    /** 
     * METHOD: operator==()
     * 
     * Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {return false;}
      if (this->node1_ == e.node1_ && this->node2_ == e.node2_) {return true;}
      if (this->node1_ == e.node2_ && this->node2_ == e.node1_) {return true;}
      return false;
    }



    /** 
     * METHOD: operator<()
     * 
     * Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // If in same graph, compare the endpoint nodes' UIDs lexicographically.
      if (this->graph_==e.graph_) {
        if (this->node1_ < this->node2_) {
          if (e.node1_ < e.node2_) {
            if (this->node1_ < e.node1_) {return true;}
            else if(this->node1_ == e.node1_) {return this->node2_ < e.node2_;}
            else {return false;}
          }
          else {
            if (this->node1_ < e.node2_) {return true;}
            else if(this->node1_ == e.node2_) {return this->node2_ < e.node1_;}
            else {return false;}
          }
        }
        else {
          if (e.node1_ < e.node2_) {
            if (this->node2_ < e.node1_) {return true;}
            else if(this->node2_ == e.node1_) {return this->node1_ < e.node2_;}
            else {return false;}
          }
          else {
            if (this->node2_ < e.node2_) {return true;}
            else if(this->node2_ == e.node2_) {return this->node1_ < e.node1_;}
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
    friend class Graph::IncidentIterator;
    friend class Graph::EdgeIterator;

    Edge(const graph_type* graph_pointer,
         size_type node1_uid,
         size_type node2_uid)
            : graph_(graph_pointer),
              node1_(node1_uid),
              node2_(node2_uid) {}

    // Pointer to parent graph
    const graph_type* graph_;

    // Node1's UID
    size_type node1_;

    // Node2's UID
    size_type node2_;
  };



  /** 
   * METHOD: num_edges()
   * 
   * Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // count edges seen until now
    size_type count = 0;

    for(auto iter = edge_map_ord_values_.begin();
        iter != edge_map_ord_values_.end();
        ++iter) {
      count += (iter->second).size();
    }
    return count;
  }


  /** 
   * METHOD: edge()
   * 
   * Return the edge with index @a i,
   * or invalid node if i >= num_edges()
   * 
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  Edge edge(size_type i) const {
    // count edges seen until now
    size_type count = 0;

    for(auto iter = edge_map_ord_values_.begin();
        iter != edge_map_ord_values_.end();
        ++iter) {
      if (count + (iter->second).size() > i) {
        return Edge(this,
                    iter->first,
                    (iter->second)[i-count].first);
      }
      count += (iter->second).size();
    }

    // We only reach this line if i >= num_edges()
    assert(false);
    return Edge();
  }



  /** 
   * METHOD: has_edge()
   * 
   * Test whether two nodes are connected by an edge.
   * 
   * @pre @a a and @a b are valid nodes of this graph
   * 
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i=0; i<edge_map_unord_.at(a.uid_).size(); i++) {
      if (edge_map_unord_.at(a.uid_)[i] == b.uid_) { return true; }
    }
    return false;
  }



  /** 
   * METHOD: add_edge()
   * 
   * Add an edge to the graph, or return the current edge if it already exists.
   * 
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * 
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()).
   */
  Edge add_edge(const Node& a, const Node& b) {
    Edge new_edge(this, a.uid_, b.uid_);
    for (size_type i=0; i<edge_map_unord_.at(a.uid_).size(); i++) {
      if (edge_map_unord_.at(a.uid_)[i] == b.uid_) return new_edge;
    }

    //Edge not present in graph
    edge_map_unord_.at(a.uid_).push_back(b.uid_);
    edge_map_unord_.at(b.uid_).push_back(a.uid_);
    if (a.uid_ <= b.uid_) {
      edge_map_ord_values_.at(a.uid_).push_back(
              std::make_pair(b.uid_, edge_value_type())
              );
    } else {
      edge_map_ord_values_.at(b.uid_).push_back(
              std::make_pair(a.uid_, edge_value_type())
              );
    }

    return new_edge;
  }





  /** 
   * METHOD: remove_edge() [given its endpoints]
   * 
   * @param[in] n1, node in the graph
   * @param[in] n2, node in the graph
   * 
   * Remove edge (n1,n2) from its graph.
   * 
   * @pre n1 is a valid node in its graph
   * @pre n2 is a valid node in its graph
   * @pre n1 and n2 are in the same graph
   * 
   * @post if (n1,n2) was not an edge, graph is not modified
   *       and new num_edges() == old num_edges()
   * @post if (n1,n2) was an edge, it's been removed
   *       and new num_edges() == old num_edges() - 1
   * @post if Edge e has e.node1()==n1 and e.node2()==n2,
   *       or viceversa, then ! new has_edge(e)
   * 
   * INVALIDATION SPECIFICATIONS:
   * @post Doesn't invalidate any node objects.
   * @post Invalidates an edge object e only if e ==
   *       the node we're removing.
   * @post Doesn't invalidate any nodeiterators.
   * @post Doesn't invalidate any incidentiterators unless
   *       their source node is one of the endpoints of the
   *       edge we're removing.
   * @post Doesn't invalidate edgeiterators unless their first
   *       node is one of the endpoints of the edge we're
   *       removing.
   * 
   * @return 0 if no edge removed, else 1
   * 
   * Complexity: O(degree(n1) + degree(n2)).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // If nodes not in graph, do nothing
    if (!has_node(n1) || !has_node(n2)) { return 0; }

    // Remove from edge_map_unord_ n1 entry
    bool removed = false;
    for (auto it = edge_map_unord_.at(n1.uid_).begin();
         it != edge_map_unord_.at(n1.uid_).end();
         ++it) {
      if (*it == n2.uid_) {
        removed = true;
        edge_map_unord_.at(n1.uid_).erase(it);
        break;
      }
    }
    // If not removed, edge wasn't in graph.
    // No need to modify anything
    if (!removed) { return 0; }

    // Remove from edge_map_unord_ n2 entry
    removed = false;
    for (auto it = edge_map_unord_.at(n2.uid_).begin();
         it != edge_map_unord_.at(n2.uid_).end();
         ++it) {
      if (*it == n1.uid_) {
        removed = true;
        edge_map_unord_.at(n2.uid_).erase(it);
        break;
      }
    }
    assert(removed);

    // Remove from edge_map_ord_values_
    size_type smaller = (n1.uid_ <= n2.uid_) ? n1.uid_ : n2.uid_;
    size_type larger = (n1.uid_ <= n2.uid_) ? n2.uid_ : n1.uid_;

    removed = false;
    for (auto it = edge_map_ord_values_.at(smaller).begin();
         it != edge_map_ord_values_.at(smaller).end();
         ++it) {
      if (it->first == larger) {
        removed = true;
        edge_map_ord_values_.at(smaller).erase(it);
        break;
      }
    }
    assert(removed);
    return 1;
  }


  /** 
   * METHOD: remove_edge() [given the edge itself]
   * 
   * Remove e from its graph, if it is still there.
   * 
   * @param[in] e, an edge. 
   * 
   * @pre e's endpoints are valid nodes in e's graph.
   * @pre e.graph_ points to this graph.
   * 
   * @post if e was not an edge, graph is not modified
   *       and new num_edges() == old num_edges()
   * @post if e was an edge, it's been removed
   *       and new num_edges() == old num_edges() - 1
   * @post !has_edge(e)
   * 
   * INVALIDATION SPECIFICATIONS:
   * Same as remove_edge(const Node&, const Node);
   * 
   * @return 0 if no edge removed, else 1
   * 
   * Complexity: O(degree(e.node1()) + degree(e.node2()))
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }



  /**
   * METHOD: remove_edge() [given an iterator to the edge]
   * 
   * Remove an edge, given an edgeiterator to it
   * 
   * @param[in] e_it, an edge iterator pointing to a valid
   *            edge in this graph.
   * 
   * @pre *e_it is a valid edge in this graph.
   * 
   * @post  if e_it points to a valid edge in the graph
   *        (and so, for instance, it wasn't edge_end()), then
   *        new num_edges() == old num_edges() - 1.
   *        Otherwise,
   *        new num_edges() == old num_edges()
   * 
   * INVALIDATION SPECIFICATIONS:
   * Same as remove_edge(const Node&, const Node&);
   * In particular, it invalidates e_it.
   * 
   * @return a valid edgeiterator of this graph, though
   *         it might be edge_end() if we removed the only edge.
   * 
   * Complexity: O(degree(n1)+degree(n2)),
   *             where n1 and n2 are the endpoints of *e_it
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }






  /** 
   * METHOD: clear()
   * 
   * Remove all nodes and edges from this graph.
   * Invalidates all outstanding Node and Edge objects.
   * 
   * Lose all internal info like nodes' positions,
   * edges' values, and nodes' values.
   * 
   * @post num_nodes() == 0 && num_edges() == 0
   */
  void clear() {
    nodes_.clear();
    node_uid_to_idx_.clear();
    edge_map_unord_.clear();
    edge_map_ord_values_.clear();
    /**
     * next node UID is not set to 0
     * because the graph has to remember that the nodes
     * we already removed from it (via this method) are
     * no longer members of it, even if we add new nodes
     * in the future (after calling this method)
     */
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
    using reference         = Node&;                    //Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  //Weak Category,Proxy

    /** 
     * CONSTRUCTOR: NodeIterator
     * 
     * Construct an invalid node iterator.
     * Points to no node in no graph.
     */
    NodeIterator() : graph_(nullptr), index_(0) {}

    /**
     * METHOD: operator*()
     * 
     * @pre nodeiterator points to existing graph
     * 
     * @pre nodeiterator is not pointing at the end
     *      (i.e. one past the last node),
     *      but rather it is pointing to a valid node in the graph.
     *      So, for instance, *(graph.node_end()) is illegal.
     * 
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
     * @return a reference to this same nodeiterator,
     *         after the ++ has been performed.
     */
    NodeIterator& operator++() {
      if (index_ < graph_->num_nodes()) { index_++; }
      return *this;
    }

    /**
     * METHOD: operator==()
     * 
     * @param[in] @a ni is a reference to a NodeIterator that the function
     *            does not modify
     * 
     * @return whether or not the two nodeiterators point to the same node 
     *         (or are both at the end) of the same graph.
     * 
     * @post the method does not modify *this
     */
    bool operator==(const NodeIterator& ni) const {
      return (graph_ == ni.graph_) && (index_ == ni.index_);
    }

  private:
    friend class Graph;

    // Parent graph
    const graph_type* graph_;

    // Pointed node's index in the graph
    size_type index_;

    /**
     * CONSTRUCTOR: NodeIterator
     * 
     * Given the parent graph and the node's index in the graph
     * construct a nodeiterator to that node.
     */
    NodeIterator(const graph_type* graph_pointer, size_type index) 
            : graph_(graph_pointer), index_(index) {}
  };

  /**
   * METHOD: node_begin()
   * 
   * @return a nodeiterator of the graph pointing to the first node.
   *         or to the end() in the case that the graph has zero nodes
   * 
   * @post the method does not modify *this
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
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
    return node_iterator(this, num_nodes());
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
    using iterator_category = std::forward_iterator_tag;  //Weak Category Proxy

    /** 
     * CONSTRUCTOR: IncidentIterator
     * 
     * Construct an invalid incident iterator.
     * There is no graph, no source node, and no index node.
     * an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr), source_(0), index_(0) {}

    /**
     * METHOD: operator*()
     * 
     * @pre incidentiterator belongs to an existing graph
     *                       and points to a valid edge
     *                       (i.e. is not at the end())
     * 
     * @return the Edge it points to
     * 
     * @post this method does not modify @this
     */
    Edge operator*() const {
      if (graph_==nullptr) return Edge();
      if (index_ < graph_->edge_map_unord_.at(source_).size()) {
        return Edge(graph_,
                    source_,
                    graph_->edge_map_unord_.at(source_)[index_]);
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
      if (graph_==nullptr) { return *this; }
      if (index_<graph_->edge_map_unord_.at(source_).size()) { index_++; }
      return *this;
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
     *         enough for the edges to match; the "root node" (source_)
     *         must match as well.
     *         Two invalid incidentiterators will also evaluate to being equal.
     */
    bool operator==(const IncidentIterator& ii2) const {
      return (graph_ == ii2.graph_) && 
             (source_ == ii2.source_) &&
             (index_ == ii2.index_);
    }

  private:
    friend class Graph;
    friend class Graph::Node;

    /**
     * CONSTRUCTOR: IncidentIterator
     * 
     * Construct an IncidentIterator given the graph, the source node,
     * and the index of that source node's list of adjacent edges.
     */
    IncidentIterator(const graph_type* graph_pointer,
                     size_type source,
                     size_type index) :
                     graph_(graph_pointer),
                     source_(source),
                     index_(index) {}

    // Pointer to parent graph
    const graph_type* graph_;

    // Source node's UID
    size_type source_;

    // Index over the vector of neighbours of the source node
    size_type index_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. 
   */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    //Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  //Weak Category,Proxy

    /** Construct an invalid EdgeIterator.
     * Points to no graph. the map_iterator_ is a default map iterator. 
     * The index defaults to 0
     */
    EdgeIterator() : graph_(nullptr), map_iterator_(), index_(0) {}

    /**
     * METHOD: operator*()
     * 
     * @pre The edgeiterator belongs to a valid graph
     * 
     * @post this method does not modify *this
     * 
     * @return the Edge this edgeiterator points to in its parent graph,
     *         or an invalid Edge object in the case that the edgeiterator
     *         is pointing at the end.
     *         E.G. if g is an empty graph, then *(g.edge_begin()) is an
     *         invalid Edge object.
     *         E.G. For any graph g, *(g.edge_end()) is an invalid edge object.
     */
    Edge operator*() const {
      if (map_iterator_ == graph_->edge_map_ord_values_.end()) return Edge();
      // If index is after the end of the source node's vector of neighbours,
      // return invalid edge.
      // However, this case should never happen in practice
      if (index_ >= map_iterator_->second.size()) { 
        assert(false);
        return Edge();
      }

      // Get node endpoints
      return Edge(graph_,
                  map_iterator_->first,
                  map_iterator_->second[index_].first);
    }

    /**
     * METHOD: operator++()
     * 
     * @post The edgeiterator is not modified in the case that it is an
     *       invalid edgeiterator object
     *       or in the case that it was already pointing at the end.
     *       Otherwise, the edgeiterator now points to the next available
     *       edge, or to the end in the
     *       case that it was previously pointing to the last available edge.
     * 
     * @return a reference to *this
     */
    EdgeIterator& operator++() {
      // If invalid edgeiterator object, return
      if (graph_==nullptr) { return *this; }

      // If we're at end, return 
      if (map_iterator_ == graph_->edge_map_ord_values_.end()) { return *this; }

      // If index is not at end of source node's vector of neighbours,
      // just increment index
      if (index_ < map_iterator_->second.size()-1) {
        index_++;
        return *this;
      }

      // Find next node with nonzero amount of neighbours 
      ++map_iterator_;
      while (map_iterator_ != graph_->edge_map_ord_values_.end()) {
        if (map_iterator_->second.size()==0) {
          ++map_iterator_;
          continue;
        }
        index_ = 0;
        return *this;
      }

      // If we're here, then we found no more edges. Return the end()
      map_iterator_ = graph_->edge_map_ord_values_.end();
      index_ = 0;
      return *this;
    }

    /**
     * METHOD: operator==()
     * 
     * @param[in] @a ei1, a edgeiterator object to compare this to.
     * 
     * @post this method does not modify *this
     * 
     * @return whether or not the edgeiterators point to the same edge
     *         in the same graph.
     *         Behaviour when one or both of the edgeiterator objects is
     *         invalid is undefined.
     */
    bool operator==(const EdgeIterator& ei2) const {
      return (graph_ == ei2.graph_) &&
             (map_iterator_ == ei2.map_iterator_) &&
             (index_ == ei2.index_);
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    /**
     * CONSTRUCTOR: EdgeIterator(graph_pointer)
     * 
     * @pre The graph_pointer parameter must be a valid graph.
     * 
     * Initialize an edge iterator pointing to the start given its parent graph
     */
    EdgeIterator(const graph_type* graph_pointer) : 
            graph_(graph_pointer),
            map_iterator_(graph_pointer->edge_map_ord_values_.begin()),
            index_(0) {}

    // Pointer to parent graph                                                
    const graph_type* graph_;

    // Iterator on the graph's internal edge_map_ord_values_
    // map sizetype vector pair sizetype edgevaluetype iterator type.
    msvpsei_type map_iterator_;

    // Index tells us the index of the edge we're looking at
    // in the map's entry's vector.
    size_type index_;
  };

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
    while (ei.map_iterator_ != edge_map_ord_values_.end()) {
      if (ei.map_iterator_->second.size() > 0) {
        // Found the first edge!
        ei.index_ = 0;
        return ei;
      }
      // Keep looking
      ++(ei.map_iterator_);
    }

    // If we reach here, then graph has no edges and
    // ei.map_iterator_ is at edge_map_ord_values_.end()
    ei.index_ = 0;
    return ei;
  }

  /**
   * METHOD: edge_end()
   * @pre This is a valid graph object.
   * 
   * @post the method does not modify *this
   * 
   * @return an edge iterator point to the end ("one past the last")
   *         of the graph's edges.
   *         So *(graph.edge_end()) is always an invalid Edge object.
   */
  edge_iterator edge_end() const {
    EdgeIterator ei(this);
    ei.map_iterator_ = edge_map_ord_values_.end();
    ei.index_ = 0;
    return ei;
  }

private:

  /**
   * INTERNAL STRUCTURE: internal_node
   * 
   * Store's a node UID, index, position, and value.
   * Proxy nodes read the index, position, and value from
   * this struct. However, they don't read the UID from
   * this struct, since they themselves store the UID.
   */

  struct internal_node {
    size_type uid_;
    size_type idx_;
    Point position_;
    node_value_type val_;

    internal_node(size_type uid,
                  Point position,
                  size_type idx,
                  node_value_type val)
      : uid_(uid), idx_(idx), position_(position), val_(val) {}

  };

  /**
   * ATTRIBUTE: nodes_
   * 
   * Vector of internal nodes. Stores uid, index, position, and value info.
   * Index is the actual index of the node in the vector, but it's stored
   * in the internal node structure so that proxy nodes also know the 
   * node's index in O(1) time without having to traverse this vector.
   */
  std::vector<internal_node> 
          nodes_;

  /**
   * ATTRIBUTE: next_node_uid_
   * 
   * Next available node uid. These are not recycled even if we remove nodes.
   * Thus, this variable is always strictly increasing when we add nodes.
   */
  size_type
          next_node_uid_;

  /**
   * ATTRIBUTE: node_uid_to_idx
   * 
   * Map a node's UID to its index in the graph.
   * If the node has been removed, its entry is deleted
   * from this map.
   */
  std::unordered_map<size_type, size_type>
          node_uid_to_idx_;

  /**
   * ATTRIBUTE: edge_map_unord_
   * 
   * Map a node's UID to a vector of the UIDs of its neighbours.
   * If a node has been removed, its entry is deleted from this map,
   * and for each old neighbour it had, their entries in this map
   * also remove the node from the vector of neighbours.
   * 
   * A neighbourless node maps to an empty vector.
   */
  std::map<size_type, std::vector<size_type>>
          edge_map_unord_;

  /**
   * ATTRIBUTE: edge_map_ord_values_
   * 
   * Map a node's UID to a vector of (UID, E) pairs where E is the
   * edge value of the edge (node1, node2), and node2 is a neighbour
   * of node1 SUCH THAT node1<=node2. A neighbourless node, or one
   * whose neighbours all have UID's smaller than it, maps to an
   * empty vector.
   * 
   * Thus, each node's entry in this map might not contain full info
   * about its neighbours; for that, we need to look at the unord_ map.
   * 
   * If a node is deleted, its entry in this map is removed from this map,
   * and for each old neighbour it had, their entries in this map
   * also remove the ndoe from the vector of neighbours.
   */
  std::map<size_type, std::vector<std::pair<size_type, edge_value_type>>>
            edge_map_ord_values_; 
};

#endif // CME212_GRAPH_HPP
