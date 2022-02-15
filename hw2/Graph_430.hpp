#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */


#define MAX_NODES_INT 92682 // if size_t is to be able to index the edges, this is the largest number of nodes we can contain (for a complete graph)

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <list>
#include <string>
#include <functional>
#include <limits.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
using std::vector;
using std::pair;
using std::unordered_map;
using std::list;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<class V, class E>
class Graph {
 private:

  struct internal_node;
  struct internal_edge;
  struct NodeHash;

  //
  // PUBLIC TYPE DEFINITIONS
  //
 public:

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
    // NodeIterator nitr = nodes.begin();
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
    Node () {
      parent = nullptr;
      idx = UINT_MAX;
    }

    /** Return this node's position. */
    const Point& position() const {
      return parent->nodes[idx].position;
    }

    /** Adding public function to make node position modifiable */
    Point& position() 
    {
      return parent->nodes[idx].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    /** Node's value data member */
    node_value_type& value() {
      internal_node& node = parent->nodes[idx];
      return node.value;
    }

    /** Get the Node's value data member */
    const node_value_type& value() const {
      const internal_node& node = parent->nodes[idx];
      return node.value;
    }

    /** Give the degree of a Node */
    size_type degree() const {
      return parent->incidents.at(*this).size();
    }

    /** Define start of incident iterator (i.e. for a given node, how do we define start) */
    incident_iterator edge_begin() const {
      IncidentIterator begin = IncidentIterator(parent, *this, parent->incidents.at(*this).begin());
      return begin;
    }

    /** Define end of incident iterator (i.e. for a given node, how do we define end) */
    incident_iterator edge_end() const {
      IncidentIterator end = IncidentIterator(parent, *this, parent->incidents.at(*this).end());
      return end;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (parent == n.parent && idx==n.idx);
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
      return idx < n.idx;
    }

    const Graph* get_parent() const {
      return parent;
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Declaring private data members and methods for Node
    size_type idx;
    Graph* parent;
    Node(size_type idx_,Graph* parent_): idx(idx_), parent(parent_) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodes.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    if (nodes.size() == MAX_NODES_INT) // verifying sufficient storage to add new node
      return Node();
    Node new_node = Node(nodes.size(),this);
    internal_node new_internal; // adding new internal node to struct
    incidents[new_node]; // adding all nodes to the incidents map so degree 0 nodes are included in degree function
    new_internal.value = v; // defining internal_node struct attributes
    new_internal.position = position;
    new_internal.node = new_node;
    nodes.push_back(new_internal); // adding internal node to vector of nodes
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.parent && n.index() < nodes.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return nodes[i].node;
  }
  
  /** Method to remove a node from the graph */
  size_type remove_node(const Node& n1) {
    // check if graph has node
    if (!(has_node(n1))) {
      return false;
    }
    // remove edge wherever n1 is incident
    for (NodeIterator n_it = node_begin(); n_it != node_end(); n_it++) {
      // use remove edge function to remove other occurences from incidents map and edges
      if (has_edge(n1,*n_it)) {
        remove_edge(n1,*n_it);
        remove_edge(*n_it,n1);
        continue;
      }
      else if (has_edge(*n_it,n1)) { // should already be accounted for in first if but double checking
        remove_edge(*n_it,n1);
        continue;
      }
    }
    size_type updated_idx = n1.idx; // finding index to be replaced
    std::swap(nodes[n1.idx],nodes[nodes.size()-1]); // swapping node to be removed to end
    nodes.pop_back(); // removing node in last position
    nodes[updated_idx].node.idx = updated_idx; // updating idx of previously last node to swapped idx
    return true;
  }

  /** Method to update node_iterator for node removal */
  node_iterator remove_node(node_iterator n_it ) {
    size_type updated_idx = (*n_it).idx; // get index of object attached to iterator
    std::swap(nodes[(*n_it).idx], nodes[nodes.size()-1]); // swap dereferenced n_it object to end
    nodes.pop_back(); // removing node in last position
    nodes[updated_idx].idx = updated_idx; // updating idx of previously last node to swapped idx
    return n_it; // returning valid iterator
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
      parent = nullptr;
      idx = UINT_MAX;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return parent->edges[idx].n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return parent->edges[idx].n2;
    }

    /** Edge's value data member */
    edge_value_type& value() {
      return parent->edges[idx].value;
    }

    /** Get the Edge'ss value data member */
    const edge_value_type& value() const {
      return parent->edges[idx].value;
    }

    /** function to compute euclidean distance (l2 norm) between two points */
    double length() const { 
      double distance = norm(node1().position() + (-node2().position()));
      return distance;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (parent == e.parent && e.node1()==node1() && e.node2()==node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (parent != e.parent) {
        return true;
      }
      return (node1()==e.node1() ? node2()<e.node2() : node1()<e.node1());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Declaring private data members and methods for Edge
    size_type idx;
    Graph* parent;
    Edge(size_type idx_, Graph* parent_): idx(idx_), parent(parent_) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges[i].edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (!(has_node(a) && has_node(b))) {
      return false;
    }
    else {
      return (incidents.at(a).find(b) != incidents.at(a).end());
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
  Edge add_edge(const Node& a, const Node& b) {
    if (has_edge(a,b)) {
      size_type idx = incidents[a][b];
      return edges[idx].edge; // returning existing edge if has_edge function found
    }
    else if (a.parent==this && b.parent==this) { // we want the edge to connect only vertices from the same parent graph
      Edge new_edge(edges.size(),this);
      internal_edge new_internal; // adding new internal edge struct to add to edges vector
      size_type size = num_edges(); // new "index" for edge based on current size of edges vector
      new_internal.n1 = a;
      new_internal.n2 = b;
      new_internal.edge = new_edge;
      edges.push_back(new_internal); // update edges vector for size and internal edge with private attributes
      incidents[a][b] = size; // adding edge to incidents map to be retrieved in either direction (a,b) and (b,a)
      incidents[b][a] = size;
      return new_edge;
    }
    return Edge(); // return invalid edge if conditions not met
  }

  /** Method to remove edge given two nodes */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // check if edge exists for removal, if not return false
    if (!(has_edge(n1,n2)) || !(has_edge(n2,n1)))
      return false;
    // remove from incidents map in both directions
    incidents.at(n1).erase(n2);
    incidents.at(n2).erase(n1);
    // remove from edges vector of internal_edges (swap and pop)
    for (EdgeIterator e_it = edge_begin(); e_it != edge_end(); ++e_it) {
      if (((*e_it).node1() == n1 && (*e_it).node2() == n2) 
           || ((*e_it).node1() == n2 && (*e_it).node2() == n1)) { // find edges based on nodes
        size_type updated_idx = (*e_it).idx;
        std::swap(edges[(*e_it).idx],edges[edges.size()-1]);
        edges.pop_back();
        edges[updated_idx].edge.idx = updated_idx;
        return true;
      }
    }
    return false;
  }
  
  /** Method to remove edge given edge */
  size_type remove_edge(const Edge& edge) {
    // check if edge exists for removal, if not return false
    if (!(has_edge(edge.node1(),edge.node2())) || !(has_edge(edge.node2(),edge.node1())))
      return false;
    else {
      // remove from incidents map in both directions
      incidents.at(edge.node1()).erase(edge.node2());
      incidents.at(edge.node2()).erase(edge.node1());
      // remove from edges vector (swap and pop)
      size_type updated_idx = edge.idx; 
      std::swap(edges[edge.idx],edges[edges.size()-1]);
      edges.pop_back();
      edges[updated_idx].edge.idx = updated_idx;
      return true;
    }
  }  
  
  /** Method to remove edge given edge_iterator */
  edge_iterator remove_edge (edge_iterator e_it) {
    size_type updated_idx = (*e_it).idx; // get index of object attached to iterator
    std::swap(edges[(*e_it).idx], edges[edges.size()-1]); // swap dereferenced e_it object to end
    edges.pop_back(); // removing node in last position
    edges[updated_idx].idx = updated_idx; // updating idx of previously last edge to swapped idx
    return e_it; // returning valid iterator
  }
 
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    incidents.clear();
    nodes.clear();
    edges.clear();
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
    using node_itr = typename std::vector<internal_node>::const_iterator;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    Node operator*() const {
      return (*nodes_).node; // pull the node from the dereferenced internal_node
    }

    NodeIterator& operator++() {
      nodes_++;
      return *this;
    }

    NodeIterator& operator++(int) {
      nodes_++;
      return *this;
    }
    
    bool operator==(const NodeIterator& n) const {
      return (parent_ == n.parent_ && nodes_ == n.nodes_);
    }

    bool operator!=(const NodeIterator& n) const {
      return !(parent_ == n.parent_ && nodes_ == n.nodes_);
    }

   private:
    friend class Graph;
    // declaring additional private attributes and constructor
    node_itr nodes_;
    Graph* parent_;
    NodeIterator(node_itr itr_on_graphs_nodes, Graph* parent): nodes_(itr_on_graphs_nodes), 
      parent_(parent) {}
  };

  node_iterator node_begin() const {
    NodeIterator begin = NodeIterator(nodes.begin(),const_cast<Graph*>(this));
    return begin;
  }

  node_iterator node_end() const {
    NodeIterator end = NodeIterator(nodes.end(),const_cast<Graph*>(this));
    return end;
  }

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
    using incident_itr = typename std::unordered_map<Node, size_type, NodeHash>::const_iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      size_type idx_ = (*incident_).second; // pulling edge attributes from incident_struct map
      Node node_1_ = node_;
      Node node_2_ = (*incident_).first;
      Edge new_edge = Edge(idx_, parent_);
      new_edge.node1() = node_1_;
      new_edge.node2() = node_2_;
      return new_edge; // return edge with correct ordered n1 and n2
    }

    IncidentIterator& operator++() {
      incident_++;
      return *this;
    }

    IncidentIterator& operator++(int) {
      incident_++;
      return *this;
    }
    
    bool operator==(const IncidentIterator& n) const {
      return (parent_ == n.parent_ && incident_ == n.incident_);
    }

    bool operator!=(const IncidentIterator& n) const {
      return !(parent_ == n.parent_ && incident_ == n.incident_);
    }

   private:
    friend class Graph;
    // declaring additional private attributes and constructor
    Graph* parent_;
    Node node_;
    incident_itr incident_;
    IncidentIterator(Graph* parent, Node node, incident_itr itr_on_nodes_incidents): 
      parent_(parent), node_(node), incident_(itr_on_nodes_incidents) {}
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
    using edge_itr = typename std::vector<internal_edge>::const_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    Edge operator*() const {
      return (*edge_itr_).edge;
    }
    
    EdgeIterator& operator++() {
      edge_itr_++;
      return *this;
    }

    bool operator==(const EdgeIterator& n) const {
      return (edge_itr_ == n.edge_itr_);
    }

    bool operator!=(const EdgeIterator& n) const {
      return !(edge_itr_ == n.edge_itr_);
    }

   private:
    friend class Graph;
    // declaring additional private attributes and constructor
    edge_itr edge_itr_;
    Graph* parent_;
    EdgeIterator(edge_itr edge, Graph* parent): edge_itr_(edge),parent_(parent) {} 
  };

  edge_iterator edge_begin() const {
    EdgeIterator begin = EdgeIterator(edges.begin(), const_cast<Graph*>(this));
    return begin;
  }

  edge_iterator edge_end() const {
    EdgeIterator end = EdgeIterator(edges.end(), const_cast<Graph*>(this));
    return end;
  }

 private:
  /** Internal Node structure to hold geometric position, value and actual Node object for each node */
  struct internal_node {
    Point position;
    Node node;
    node_value_type value;
    Point initial_position;
  };

  /** Internal Edge structure to hold associated nodes, edge value attribute and actual Edge object for each edge*/
  struct internal_edge {
    Node n1;
    Node n2;
    Edge edge;
    edge_value_type value;
  };

  // code for hash functions referenced from https://en.cppreference.com/w/cpp/utility/hash
  struct NodeHash
  {
      std::size_t operator()(Node const& node) const noexcept
      {
          return std::hash<size_type>{}(node.idx);
      }
  };

  /** Private data members to track and link Edges, Nodes, and their
   *  geometric positions */
  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;
  std::unordered_map<Node, std::unordered_map<Node, size_type, NodeHash>,NodeHash> incidents;
};


#endif // CME212_GRAPH_HPP