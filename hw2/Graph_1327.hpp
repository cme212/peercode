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
 */
template <typename V, typename E>
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
  /** Synonym for node value type for the graph template */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for edge value type for the graph template */
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
  Graph() {
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
  class Node: private totally_ordered <Node> {
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
      return *(graph_->node_vector[node_idx_].first);
    }

    /** Return this node's position, modifiable.*/
    Point & position () {
        return *(graph_->node_vector[node_idx_].first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(node_idx_ < graph_->size());
      return node_idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return this node's value with node_value_type. */
    node_value_type & value (){
        return graph_->node_vector[node_idx_].second;
    }

    /** Return this node's value with const node_value_type. */
    const node_value_type & value () const{
        return graph_->node_vector[node_idx_].second;
    }

    /** Return the number of nodes adjacent to the current node. */
    size_type degree() const {
        return graph_->edge_vector[node_idx_].size();
    }

    /** @return incident_iterator at the starting position */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, node_idx_, 0);
    }

    /** @return incident_iterator at one passing the ending position */
    incident_iterator edge_end() const {
        return IncidentIterator(graph_, node_idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (n.graph_ == graph_ and n.node_idx_ == node_idx_){
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
      // HW0: YOUR CODE HERE
      if (graph_ < n.graph_){
        return true;
      }
      else if (graph_ == n.graph_ and node_idx_ < n.node_idx_){
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_idx_;
    // constructor
    Node(const graph_type* graph, size_type idx):
    graph_(const_cast<graph_type*>(graph)), node_idx_(idx){}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_vector.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type & value = node_value_type()) {
    // HW0: YOUR CODE HERE
    Point* cur_pt = new Point;
    *cur_pt = position;
    std::pair<Point*, node_value_type> cur(cur_pt, value);
    node_vector.push_back(cur);
    // add new vector to the edge_vector
    std::vector<std::pair<size_type, edge_value_type>> new_vec;
    edge_vector.push_back(new_vec);
    return Node(this, size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    assert(n.graph_ != NULL);
    if(n.node_idx_ < size() and n.graph_ == this){return true;}
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
    assert(i >= 0 and i < size());
    return Node(this, i); // Valid node
  }

  /** Remove a node from the graph, returning true if sucessfully
   *  removed and false if no node to be removed.
   * @param[in] n The node to be removed
   * @post new size() == old size() - 1
   * @post new num_edges == old num_edges() - degree()
   * Complexity: O(num_nodes())
   * Invalidate: node_idx_ == n.node_idx_, the edges incident to
   * node n and the corresponding node and edge iterators.
   */
  size_type remove_node (const Node& n){
      if(!has_node(n)){
          return (size_type) false;
      }
      delete node_vector[n.node_idx_].first;
      size_type last = size() - 1;
      //std::cout << last <<"check size"<<std::endl;
      size_type aid = n.node_idx_;
      node_vector[aid] = node_vector.back();
      node_vector.pop_back();
      for(size_type i=0; i < edge_vector[aid].size(); ++i){
          size_type bid = edge_vector[aid][i].first;
          for(size_type k=0; k < edge_vector[bid].size(); ++k){
              if(edge_vector[bid][k].first == aid){
                  edge_vector[bid][k] = edge_vector[bid].back();
                  edge_vector[bid].pop_back();
              }
          }
      }
      //std::cout<<"finish removing nodes"<<std::endl;
      for(size_type i=0; i < edge_vector[last].size(); ++i){
          size_type bid = edge_vector[last][i].first;
          for(size_type k=0; k < edge_vector[bid].size(); ++k){
              if(edge_vector[bid][k].first == last){
                  edge_vector[bid][k].first = aid;
              }
          }
      }
      //std::cout<<"finish removing edges"<<std::endl;
      edge_vector[aid] = edge_vector[last];
      edge_vector.pop_back();

      return (size_type) true;
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
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type node2_idx_ = graph_->edge_vector[node1_idx_][vec2_idx_].first;
      return Node(graph_, node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      size_type node2_idx_ = graph_->edge_vector[node1_idx_][vec2_idx_].first;
      size_type node2_idx2 = e.graph_->edge_vector[node1_idx_][vec2_idx_].first;
      bool nodes = (node1_idx_ == e.node1_idx_ && node2_idx_ == node2_idx2);
      return (nodes && (graph_ == e.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      size_type node2_idx_ = graph_->edge_vector[node1_idx_][vec2_idx_].first;
      size_type node2_idx2 = e.graph_->edge_vector[node1_idx_][vec2_idx_].first;
      if (graph_ < e.graph_){
            return true;
      }else if(node1_idx_ < e.node1_idx_){
          return true;}
      else if (graph_ == e.graph_ && node1_idx_ == e.node1_idx_){
          if(node2_idx_ < node2_idx2){
              return true;
          }
      }
      return false;
    }

    /** Return a spring constant for the edge */
    double spring_const(){
        double spring_constant = graph_->edge_vector[node1_idx_].size();
        return spring_constant;
    }

    /** Return the length of the edge */
    double length() const{
        return norm(node1().position() - node2().position());
    }

    /** Return the value of the edge with edge_value_type */
    edge_value_type& value() {
        return graph_->edge_vector[node1_idx_][vec2_idx_].second;
    }

    /** Return the value of the edge with const edge_value_type */
    const edge_value_type & value () const {
        return graph_->edge_vector[node1_idx_][vec2_idx_].second;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type node1_idx_;
    size_type vec2_idx_;
    // Constructor
    Edge(const graph_type* graph, size_type idx1, size_type idx2):
            graph_(const_cast<graph_type*>(graph)), node1_idx_(idx1), vec2_idx_(idx2){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type num = 0;
    for(const auto& k : edge_vector){
        num += k.size();
    }
    size_type numEdges = num/2;
    return numEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type node_a = a.node_idx_;
    size_type node_b = b.node_idx_;
    for (size_type i = 0; i < a.degree(); i++){
      if(edge_vector[node_a][i].first == node_b){
        return true;
      }
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
    assert(a.index() != b.index());
    assert(this == a.graph_ and this == b.graph_);
    size_type id_a = a.node_idx_;
    size_type id_b = b.node_idx_;
    if (has_edge(a,b)){
        for(size_type i=0; i< a.degree(); i++){
            if(edge_vector[id_a][i].first == id_b){
                return Edge(this, id_a, i);
            }
        }
    }
    edge_value_type edgeValue = edge_value_type();
    std::pair<size_type, edge_value_type> pair_b(id_b, edgeValue);
    std::pair<size_type, edge_value_type> pair_a(id_a, edgeValue);
    edge_vector[id_a].push_back(pair_b);
    edge_vector[id_b].push_back(pair_a);
    return Edge(this, id_a, edge_vector[id_a].size()-1);
  }

  /** Remove an edge from the graph given the incident nodes.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @param[in] a,b nodes incident to the edge
   * @return size_type true if successfully removed, false otherwise
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges()-1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate the edge that connects node @a a and node @b b.
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge (const Node & a, const Node & b){
      assert(has_node(a));
      assert(has_node(b));
      if(!has_edge(a, b)){
          return (size_type) false;
      }
      size_type id_a = a.node_idx_;
      size_type id_b = b.node_idx_;
      // remove from the first end
      for(size_type i=0; i < edge_vector[id_a].size(); i++){
          if(edge_vector[id_a][i].first == id_b){
              size_type last = edge_vector[id_a].size()-1;
              edge_vector[id_a][i] = edge_vector[id_a][last];
              edge_vector[id_a].pop_back();
          }
      }
      // remove from the other end
      for(size_type i=0; i < edge_vector[id_b].size(); i++){
          if(edge_vector[id_b][i].first == id_a){
              size_type last = edge_vector[id_b].size()-1;
              edge_vector[id_b][i] = edge_vector[id_b][last];
              edge_vector[id_b].pop_back();
          }
      }
      return (size_type) true;
  }

  /** Remove an edge from the graph given the edge.
   * @param[in] e The edge to be removed
   * @return size_type true if successfully removed, false otherwise
   * @post If old has_edge(e.node1(), e.node2()),
   *            new num_edges() == old num_edges()-1.
   *       Else,
   *            new num_edges() == old num_edges().
   *
   * Can invalidate the the given edge @e e.
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge (const Edge & e){
      return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph given the edge.
   * @param[in] e_it The edge_iterator that points to the edge
   * @return the edge_iterator after the edge being removed
   * @post If old has_edge(e.node1(), e.node2()),
   *            new num_edges() == old num_edges()-1.
   *       Else,
   *            new num_edges() == old num_edges().
   *
   * Can invalidate the the edge pointed by the given edge_iterator.
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge (edge_iterator e_it){
      remove_edge(*e_it);
      e_it.reset();
      return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // delete the points allocated on heap
    for(auto i=node_vector.begin(); i!=node_vector.end(); ++i){
        delete i->first;
    }
    node_vector.clear();
    edge_vector.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered <NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag; // Weak Category, Proxy
    //using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /**
     * @brief Dereference operator for Node Iterator
     * @return Node at the position of corresponding Node Iterator
    */
    Node operator*() const {
        return Node(graph_, node_idx_);
    }

    /**
     * @brief Increments the Node Iterator to the next position
     * @return NodeIterator with position advanced by one
    */
    NodeIterator& operator++() {
        node_idx_++;
        return *this;
    }

    /**
     * @brief Test if two Node Iterators are equal
     * @return true when both Node Iterators point to the same graph
     * and at the same position
    */
    bool operator==(const NodeIterator& other) const {
        return (other.graph_ == graph_ && other.node_idx_ == node_idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type node_idx_;
    // constructor
    NodeIterator(const graph_type* graph, size_type idx):
              graph_(const_cast<graph_type*>(graph)), node_idx_(idx){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** @return node_iterator at the starting position */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** @return node_iterator at one passing the ending position */
  node_iterator node_end() const {
      return NodeIterator(this, num_nodes());
  }

  /** Remove the pointed node given the node_iterator.
   * @param[in] n_it The node_iterator
   * @post new graph_.size() == old graph_.size() - 1
   * @post new graph_.num_edges == old graph_.num_edges() - graph_.degree()
   * Complexity: O(num_nodes())
   * Invalidate: The given node_iterator.
   */
  node_iterator remove_node (node_iterator n_it){
      remove_node(*n_it);
      return n_it;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator> {
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

    /**
     * @brief Dereference operator for Incident Iterator
     * @return Edge at the position of corresponding Incident Iterator
     */
    Edge operator*() const {
        return Edge(graph_, node1_idx_, vec2_idx_);
    }

    /**
     * @brief Increments the Incident Iterator to the next position
     * @return IncidentIterator with position advanced by one
     */
    IncidentIterator& operator++() {
        vec2_idx_++;
        return *this;
    }

    /**
     * @brief Test if two Incident Iterators are equal
     * @return true when both Incident Iterators point to the same graph
     * and at the same position and same edge connecting the same nodes
     */
    bool operator==(const IncidentIterator& other) const {
        return (graph_ == other.graph_ && node1_idx_ == other.node1_idx_
            && vec2_idx_ == other.vec2_idx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node1_idx_;
    size_type vec2_idx_;

    // Constructor
    IncidentIterator(const graph_type* graph, size_type idx1, size_type idx2):
              graph_(const_cast<graph_type*>(graph)), node1_idx_(idx1), vec2_idx_(idx2){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered <EdgeIterator> {
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

    /**
     * @brief Dereference operator for Edge Iterator
     * @return Edge at the position of corresponding Edge Iterator
     */
    Edge operator*() const {
        return Edge(graph_, node1_idx_, vec2_idx_);
    }

    /**
     * @brief Increments the Edge Iterator to the next position
     * @return EdgeIterator with position advanced by one
     */
    EdgeIterator& operator++() {
        vec2_idx_++;
        size_type size1 = graph_->edge_vector.size();
        while(node1_idx_ < size1){
            size_type size2 = graph_->edge_vector[node1_idx_].size();
            while(vec2_idx_ < size2){
                size_type idx2 = graph_->edge_vector[node1_idx_][vec2_idx_].first;
                if(node1_idx_ < idx2){
                    return *this;
                }
                ++vec2_idx_;
            }
            vec2_idx_= 0;
            ++node1_idx_;
        }
        return *this;
    }

    /**
     * @brief Test if two Edge Iterators are equal
     * @return true when both Edge Iterators point to the same graph
     * and at the same node position and vector position.
     */
    bool operator==(const EdgeIterator& other) const {
        return (graph_ == other.graph_ && node1_idx_ == other.node1_idx_
                && vec2_idx_ == other.vec2_idx_);
    }

    /**
     * @brief Helper function to reset the index of the edge_iterator
     */
    void reset(){
        size_type size1 = graph_->edge_vector.size();
        while(node1_idx_ < size1){
            size_type size2 = graph_->edge_vector[node1_idx_].size();
            while(vec2_idx_ < size2){
                size_type idx2 = graph_->edge_vector[node1_idx_][vec2_idx_].first;
                if(node1_idx_ < idx2){
                    return;
                }
                ++vec2_idx_;
            }
            ++node1_idx_;
            vec2_idx_= 0;
        }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type node1_idx_;
    size_type vec2_idx_;

    // Constructor
    EdgeIterator(const graph_type* graph, size_type idx1, size_type idx2):
        graph_(const_cast<graph_type*>(graph)), node1_idx_(idx1), vec2_idx_(idx2){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** @return edge_iterator at the starting position */
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0, 0);
  }

  /** @return edge_iterator at one passing the ending position */
  edge_iterator edge_end() const {
      return EdgeIterator(this, num_nodes(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector <std::pair<Point*, node_value_type>> node_vector;
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> edge_vector;

};

#endif // CME212_GRAPH_HPP
