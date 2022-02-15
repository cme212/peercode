#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  /** Synonym for Node val type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for Edge val type. */
  using edge_value_type = E;

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
  Graph(): idToIdIndexMap_(), nodes_(), edges_() {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return *element_->point;
    }

    Point& position () {
      return *element_->point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return element_->index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return a reference of the node's value. */
    node_value_type& value() {
      return element_->val;
    }

    /** Return a const reference of the node's value. */
    const node_value_type& value() const {
      return element_->val;
    }

    /** return the number of incident edges. */
    size_type degree() const {
      auto map = element_->graph->idToIdIndexMap_.find(index());
      if (map != element_->graph->idToIdIndexMap_.end()) {
        return map->second.size();
      }
      return 0;
    }

    /** Start of the incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(element_->graph, this);
    }

    incident_iterator edge_end() const {
      IncidentIterator incident_iter = IncidentIterator(element_->graph, this);
      incident_iter.move_to_end();
      return incident_iter;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // nodes need to have the same index and belong the same graph
      if (element_->index == n.element_->index and element_->graph == n.element_->graph) {
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
     *
     * The following code compares two Point's value sequentially
     * Compare Point.x first, if the x value differs we can return immediately, then Point.y, at last Point.z
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (element_->graph < n.element_->graph) return true;
      else if (element_->graph == n.element_->graph and element_->index < n.element_->index) return true;
      else return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class EdgeIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // a struct holding essential elements
    struct internal_element {
        graph_type* graph;
        Point* point;
        size_type index;
        node_value_type val;
        internal_element(graph_type* graph, Point* point, size_type index): graph(graph), point(point), index(index) {};
        internal_element(graph_type* graph, Point* point, size_type index, node_value_type val):
        graph(graph), point(point), index(index), val(val) {};
    };
    internal_element* element_; // By using this pointer to a struct, we save some space

    Node(graph_type* graph, Point* point, size_type index) {element_ = new internal_element(graph, point, index); }
    Node(graph_type* graph, Point* point, size_type index, node_value_type val) {
      element_ = new internal_element(graph, point, index, val);
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
//  Node add_node(const Point& position) {
//    // HW0: YOUR CODE HERE
//    auto* point = new Point;
//    *point = position;
//    Node node = Node(this, point, num_nodes());
//    nodes_.push_back(node);
//    return node;
//  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type()) {
    auto* point = new Point;
    *point = position;
    Node node = Node(this, point, num_nodes(), node_val);
    nodes_.push_back(node);
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.element_->graph == this and n.element_->index < size()) return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return nodes_.at(i);
  }

  /** Remove a given node
   * Return the number of nodes remove (0, 1)
   */
   size_type remove_node(const Node& n) {
     if (!has_node(n)) {
       return 0;
     }

     // Remove the node's edges
     std::vector<edge_type> edges;
     for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
       edges.push_back(*it);
     }

     for (auto& edge : edges) {
       remove_edge(edge);
     }

     // Remove the actual node
     size_type old_size = size();
     size_type old_index = n.index();
     size_type old_degree = nodes_.back().degree();
     nodes_.at(old_index) = nodes_.back();
     nodes_.at(old_index).element_->index = old_index;
     nodes_.pop_back();

     // Replace the previous end node's index
     if (num_nodes() > 0 and old_index != old_size - 1 and old_degree != 0) {
       node_type node = nodes_.at(old_index);
       for (auto it = node.edge_begin(); it != node.edge_end(); ++it) {
         edge_type edge = *it;
         if (edge.elements_->node1_idx_ == old_size - 1) {
           edge.elements_->node1_idx_ = old_index;
         } else {
           edge.elements_->node2_idx_ = old_index;
         }
       }

       idToIdIndexMap_[old_index] = idToIdIndexMap_.at(old_size - 1);
       idToIdIndexMap_.erase(old_size - 1);

       for (auto it = idToIdIndexMap_.at(old_index).begin(); it != idToIdIndexMap_.at(old_index).end(); ++it) {
         auto other_idx = (*it).first;
         auto edge_idx = (*it).second;
         idToIdIndexMap_.at(other_idx).erase(old_size - 1);
         idToIdIndexMap_.at(other_idx).insert({old_index, edge_idx});
         edgeIndexMap_.at(edge_idx) = {old_index, other_idx};
       }
     }

     return 1;
   }

  /** Remove a given node using a node iterator
  * Return a node iterator equals to node_begin()
  */
  size_type remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return elements_->graph_->node(elements_->node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return elements_->graph_->node(elements_->node2_idx_);
    }

    /** Return the length of the edge */
    double length() const {
      return norm(elements_->bigNode_->position() - elements_->smallNode_->position());
    }

    /** Return the reference of the Edge's value */
    edge_value_type& value() {
      return elements_->value_;
    }

    const edge_value_type& value() const {
      return elements_->value_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (elements_->graph_ == e.elements_->graph_ and
         *(elements_->smallNode_) == *(e.elements_->smallNode_) and
         *(elements_->bigNode_) == *(e.elements_->bigNode_)) return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (elements_->graph_ < e.elements_->graph_ or
          elements_->smallNode_ < e.elements_->smallNode_ or
          (elements_->smallNode_ == elements_->smallNode_ and elements_->bigNode_ < elements_->bigNode_)) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    struct edge_element {
        const graph_type* graph_;
        Node* smallNode_;
        Node* bigNode_;
        size_t node1_idx_;
        size_t node2_idx_;
        edge_value_type value_;

        edge_element(const graph_type* graph, Node& a, Node& b): node1_idx_(a.index()), node2_idx_(b.index()) {
          // nodes must be on the same graph and two nodes cannot be the same
          assert(a.element_->graph == b.element_->graph and a.element_->index != b.element_->index);
          a = graph->node(a.index());
          b = graph->node(b.index());
          if (a < b) {smallNode_ = &a; bigNode_ = &b;}
          else {smallNode_ = &b; bigNode_ = &a;}
          graph_ = graph;
        }

//        edge_element(const graph_type* graph, Node* a, Node* b) {
//          // nodes must be on the same graph and two nodes cannot be the same
//          assert(a->element_->graph == b->element_->graph and a->element_->index != b->element_->index);
//          a = graph->node(a->index());
//          b = graph->node(b->index());
//          if (a < b) {smallNode_ = a; bigNode_ = b;}
//          else {smallNode_ = b; bigNode_ = a;}
//          graph_ = graph;
//        }

    };

    edge_element* elements_;

    Edge(const graph_type* graph, Node& a, Node& b) {
      elements_ = new edge_element(graph, a, b);
      set_node1_node2(a, b);
    }

//    Edge(const graph_type* graph, Node* a, Node* b) {
//      elements_ = new edge_element(graph, a, b);
//      set_node1_node2(a, b);
//    }

    /** set the edge's node1 and node2 */
    void set_node1_node2(const Node& node1, const Node& node2) {
        elements_->node1_idx_ = node1.index();
        elements_->node2_idx_ = node2.index();
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return edges_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (idToIdIndexMap_.find(a.index()) != idToIdIndexMap_.end()) {
        auto connectedNodes = idToIdIndexMap_.find(a.index())->second;
        if (connectedNodes.find(b.index()) != connectedNodes.end()) return true;
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
    if (has_edge(a, b)) {
        Edge cur_edge = edge(idToIdIndexMap_.find(a.index())->second.find(b.index())->second);
        cur_edge.set_node1_node2(a, b);
        return cur_edge;
    }
    idToIdIndexMap_[a.index()][b.index()] = edges_.size();
    idToIdIndexMap_[b.index()][a.index()] = edges_.size();
    edgeIndexMap_[edges_.size()] = {a.index(), b.index()};
    Edge edge = Edge(this, const_cast<Node&>(a), const_cast<Node&>(b));
    edges_.push_back(edge);
    return edge;
  }

  /** Remove the edge connecting two given nodes
   * Return the number of edges removed (0, 1)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) {
      return 0;
    }

    // remove the actual edge
    size_type edge_idx = idToIdIndexMap_.at(a.index()).at(b.index());
    edges_.at(edge_idx) = edges_.back();
    edges_.pop_back();

    // remove their index
    remove_edge_idx(a.index(), b.index());
    remove_edge_idx(b.index(), a.index());

    // adjust previous end edge's index
    edgeIndexMap_.erase(edge_idx);
    if (!edges_.empty()) {
      size_type node1_idx = edgeIndexMap_[edges_.size()].first;
      size_type node2_idx = edgeIndexMap_[edges_.size()].second;
      idToIdIndexMap_[node1_idx][node2_idx] = edge_idx;
      idToIdIndexMap_[node2_idx][node1_idx] = edge_idx;

      edgeIndexMap_.erase(edges_.size());
      edgeIndexMap_[edge_idx] = {node1_idx, node2_idx};
    }

    return 1;
  }

  /** Remove the given edge
   * Return the number of edges removed (0, 1)
   */
  size_type remove_edge(const Edge& edge) {
    return remove_edge(edge.node1(), edge.node2());
  }

  /** Remove the edge given an edge iterator
   * Return an edge iterator equal to edge_begin()
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }

  void remove_edge_idx(size_type center_node_idx, size_t other_node_idx) {
    idToIdIndexMap_.at(center_node_idx).erase(other_node_idx);
    if (idToIdIndexMap_.at(center_node_idx).empty()) {
      idToIdIndexMap_.erase(center_node_idx);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    idToIdIndexMap_.clear();
    nodes_.clear();
    edges_.clear();
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const {
      return *nodes_iter_;
    }

    NodeIterator& operator++() {
      ++nodes_iter_;
      return *this;
    }

    NodeIterator& operator--() {
      --nodes_iter_;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
      if (graph_ == node_iter.graph_ and nodes_ == node_iter.nodes_ and nodes_iter_ == node_iter.nodes_iter_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const graph_type* graph_;
    std::vector<Node>* nodes_;
    typename std::vector<Node>::iterator nodes_iter_;
    explicit NodeIterator(const graph_type* graph): graph_(graph), nodes_(const_cast<std::vector<Node>*>(&graph->nodes_)),
    nodes_iter_(nodes_->begin()) {}

    /** move the iterator one pass the end */
    void move_to_end() {
      nodes_iter_ = nodes_->end();
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
    return NodeIterator(this);
  }

  node_iterator node_end() const {
    node_iterator node_iter = NodeIterator(this);
    node_iter.move_to_end();
    return node_iter;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() = default;

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
      size_t edge_idx = incident_itr->second;
      Edge edge = graph_->edges_.at(edge_idx);
      size_t node2_idx = incident_itr->first;
      auto& center_ref = const_cast<node_type&>(graph_->nodes_.at(center_idx_));
      auto& node_ref = const_cast<node_type&>(graph_->nodes_.at(node2_idx));
      edge.set_node1_node2(center_ref, node_ref);
      return edge;
    }

    IncidentIterator& operator++() {
      ++incident_itr;
      return *this;
    }

    bool operator==(const IncidentIterator& incident_iter) const {
      if (degrees_ == 0 and incident_iter.degrees_ == 0) {
        return (graph_ == incident_iter.graph_ and center_ == incident_iter.center_);
      } else if (degrees_ == 0 or degrees_ == 0) {
        return false;
      }

      if (graph_ == incident_iter.graph_ and
          *center_ == *(incident_iter.center_) and
          incident_itr == incident_iter.incident_itr) {
        return true;
      }
      return false;
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const graph_type* graph_;
    const node_type* center_;
    const int center_idx_;
    std::unordered_map<size_type, size_type>::iterator incident_itr;
    std::unordered_map<size_type, size_type>::iterator incident_end;
    size_t degrees_;
    IncidentIterator(const graph_type* graph, const node_type* node): graph_(graph), center_(node), center_idx_(node->index()) {
      auto edge_map = const_cast<graph_type*>(graph_)->idToIdIndexMap_.find(center_idx_);
      // if the node has edges
      if (edge_map != graph_->idToIdIndexMap_.end()) {
        incident_itr = edge_map->second.begin();
        incident_end = edge_map->second.end();
        degrees_ = edge_map->second.size();
      } else {
        // if the node has no edge
        degrees_ = 0;
      }
    }

    /** move the iterator one pass the end */
    void move_to_end() {
      incident_itr = incident_end;
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
    Edge operator*() const {
      return *edge_itr;
    }

    EdgeIterator& operator++() {
      ++edge_itr;
      return *this;
    }

    bool operator==(const EdgeIterator& edge_iter) const {
      return graph_ == edge_iter.graph_ and edge_itr == edge_iter.edge_itr;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const graph_type* graph_;
    typename std::vector<Edge>::iterator edge_itr;
    explicit EdgeIterator(const graph_type* graph): graph_(graph), edge_itr(const_cast<graph_type*>(graph)->edges_.begin()) {}

    void move_to_end() {
      edge_itr = const_cast<graph_type*>(graph_)->edges_.end();
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return EdgeIterator(this);
  }

  edge_iterator edge_end() const {
    edge_iterator edge_iter = EdgeIterator(this);
    edge_iter.move_to_end();
    return edge_iter;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // {node1_idx, {node2_idx, edge_idx}}
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> idToIdIndexMap_;
  // {edge_idx, {node1_idx, node2_idx}}
  std::unordered_map<size_type, std::pair<size_type, size_type>> edgeIndexMap_;
  std::vector<Node> nodes_;
  std::vector<Edge> edges_;
};

#endif // CME212_GRAPH_HPP
