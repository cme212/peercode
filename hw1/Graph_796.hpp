#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <utility> 
#include <iterator>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
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

  /** Decleration of node value type, i.e. what it represents */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    //No need to put anything in here, nodes gets automoatically instantiated
    //just as a vector with no elements.
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
      //Nothing to do here, only valid constructor should be private
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return this->graph->nodes[this->uid].point;
    }

    /** Allows Nodes value to be read and modified */
    node_value_type& value(){
      return this->graph->nodes[this->uid].value;
    }

    /** Allows Node's value to be read */
    const node_value_type& value() const {
      return this->graph->nodes[this->uid].value;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->uid;
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
      // HW0: YOUR CODE HERE
      return (this->graph == n.graph) && (this->uid == n.uid);
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
      if(this->graph == n.graph){
          if(this->uid < n.uid){
              return true;
          } else {
              return false;
          }
      } else {
          std::vector<graph_type*> pointers = {this->graph, n.graph};
          std::sort(pointers.begin(), pointers.end());
          if(pointers[0] == this->graph){
              return true;
          }
      }
      return false;
    }

    /** Returns the number of incident edges */
    size_type degree() const {
        return this->graph->nodes[this->uid].connections.size();
    }
  
    /** Start of the incident iterator for this specific node
     * @return Incident_Iterator pointing to beginning of connections vector
    */
    incident_iterator edge_begin() const {
        return IncidentIterator(this->graph->nodes[this->uid].connections.begin(), 
 				*this, this->graph);
    }
   
    /** End of the incident iterator
     * @return Incident_iterator pointing to one past end of connections vector
    */
    incident_iterator edge_end() const {
        return IncidentIterator(this->graph->nodes[this->uid].connections.end(),
				*this, this->graph);
    }
    

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph;
    size_type uid;

    Node(const graph_type* graph_, size_type uid_)
        : graph(const_cast<graph_type*>(graph_)), uid(uid_) {
    }
      

  };  


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->nodes.size();
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
  Node add_node(const Point& position,
		const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    internal_node i_node;
    i_node.uid = (size_type) this->nodes.size();
    i_node.point = position;
    i_node.value = value;

    nodes[i_node.uid] = i_node;
    node_keys.insert(i_node.uid);
    return Node(this, this->nodes.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    if(this->node_keys.find(n.uid) != this->node_keys.end()){
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
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    //We'll return a proxy object for element @i a:
    assert(i < nodes.size());
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return this->graph->node(a_index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return this->graph->node(b_index);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return (this->node1() == e.node1()) && (this->node2() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //Just piggyback off ordering of nodes
      (void) uid; //silence compiler warning, but is prob redundant
      if(this->node1() == e.node1()){
          return (this->node2() < e.node2());
      }
      return (this->node1() < e.node1());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(size_type a_index_, size_type b_index_, size_type uid_, \
         const graph_type* graph_)
        : a_index(a_index_), b_index(b_index_), uid(uid_), \
               graph(const_cast<graph_type*>(graph_)) {
    }
 
    size_type a_index;
    size_type b_index;
    size_type uid;
    graph_type* graph; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
 
    assert(i < edges.size());
    std::pair<size_type, size_type> key = edge_keys.at(i);
    return Edge(key.first, key.second, i, this);       
  }

  /** Return the edge with node index @a a and node index @a ab
   * @pre (@a a, @a b) elements of edge_keys
   * Complexity: Roughly O(1)
   * @return Edge with node1 index @a a and node2 index @a b
  */
  Edge edge(size_type a_uid, size_type b_uid) const {
    std::pair<size_type, size_type> key1 = {a_uid, b_uid};
    std::pair<size_type, size_type> key2 = {b_uid, a_uid};
    //We have to find the "correct" (a,b) key. Implementation only stores
    //(a,b) key, so querying (b,a) will return false statement even if undirected graph
    internal_edge old_edge;
    if(edges.count(key1)){
        old_edge = edges.at(key1);
	return Edge(a_uid, b_uid, old_edge.uid, this);
    }
    //We still need to return the incident node in the as the b node
    old_edge = edges.at(key2);
    return Edge(a_uid, b_uid, old_edge.uid, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    std::pair<size_type, size_type> key1 = {a.index(), b.index()};
    std::pair<size_type, size_type> key2 = {b.index(), a.index()};
    int b1 = edges.count(key1);
    int b2 = edges.count(key2);
    if (b1 == 0 && b2 ==0){
        return false;
    } 
    return true;
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
    //(void) a, (void) b;   // Quiet compiler warning
    if(has_edge(a, b)){
        std::pair<size_type, size_type> key1 = {a.index(), b.index()};
        std::pair<size_type, size_type> key2 = {b.index(), a.index()};
        internal_edge old_edge;
        if(edges.count(key1)){
            old_edge = edges.at(key1);
            return Edge(a.index(), b.index(), old_edge.uid, this);
        }
        //need to recover the id essentially via the internal node
        old_edge = edges.at(key2);
        return Edge(a.index(), b.index(), old_edge.uid, this);
    } 
    
    std::pair<size_type, size_type> new_key = {a.index(), b.index()};
    edge_keys[this->edges.size()] = new_key;
    internal_edge edge;
    edge.uid = (size_type) this->edges.size();
    edge.a = a;
    edge.b = b;
    edges[new_key] = edge;

    //add connections to internal_node
    nodes[a.index()].connections.push_back(b.index());
    nodes[b.index()].connections.push_back(a.index());

    return Edge(a.index(), b.index(), this->edges.size()-1, this);   
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    node_keys.clear();
    edges.clear();
    edge_keys.clear();

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {};

    /** Iterate to next node in nodes vector.
      * @return Reference to node iterator
    */
    NodeIterator& operator++() {
        node_itr++;
	return *this;
    }

    /** Equality for NodeIterator class 
     * @return bool, whether or not NodeIter refers to same graph and vector<Node>::iterator
     */
    bool operator==(const NodeIterator& other_iter) const{
        return (node_itr == other_iter.node_itr) &&
		 (graph == other_iter.graph);
    }
    
    /** Dereference operator
     * @pre this != NodeIterator.end()
     * @return Node to which our operator points
    */
    Node operator*(){
	size_type idx = *node_itr;
	return graph->node(idx);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    //Since we implemented internal structure of nodes via map, will use a stl map iterator
    std::set<size_type>::iterator node_itr;
    const graph_type* graph;    

    //Private constructor
    NodeIterator(std::set<size_type>::iterator set_itr, 
		const graph_type* graph_) : 
			node_itr(set_itr), graph(graph_) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** begin() function for NodeIterator
   * @return NodeIterator pointing to beginning of nodes vector
  */
  NodeIterator node_begin() const {
      return NodeIterator(node_keys.begin(), this);
  //    return nodes.erase(nodes.begin(), nodes.begin());
  }

  /** end() function for NodeIterator
   * @return NodeIterator pointing to one past end of nodes vector
  */
  NodeIterator node_end() const {
      return NodeIterator(node_keys.end(), this);
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
    IncidentIterator() {
    }


    /**  Dereference operator for IncidentIterator
     * @pre this->incident_itr != incident_itr.end()
     * @return an Edge that connects root Node and another node
     * @post Edge.node1() == root_node, Edge.node2() == connected_node
     */ 
    Edge operator*() const{ 
        //First get the uid of the other outgoing node
        size_type incident_node_id = *incident_itr;
	/** Now check that edge exists, and return it via add_edge */
        if(graph->has_edge(root_node, graph->node(incident_node_id))){
            return graph->edge(root_node.uid, incident_node_id);
        } 
        /** If edge not found, bug exists in program */
        std::cerr << "Couldn't find edge that should exist between nodes" << std::endl;
        exit(1);
    }

    /** Next operator for IncidentIterator
     * @return IncidentIterator pointing to next incident edge's node
     */
    IncidentIterator& operator++(){
        incident_itr++;
        return *this;
    }

    /** Equality operator for IncidentIterator
     * @return bool, whether or not two IncidentIterators have same root,
     *         graph, and point to same connection vector element
     */
    bool operator==(const IncidentIterator& other_iter) const{
	return (incident_itr == other_iter.incident_itr) &&
	       (root_node == other_iter.root_node) && 
	        (graph == other_iter.graph);
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;

    /**Essentially, we have to iterate over the internal nodes member 
     * variable @a connections.This iterates over the incident nodes ids,
     * allowing us to access them in constant time
    */
    std::vector<size_type>::iterator incident_itr;
    Node root_node;
    const graph_type* graph;
  
    //Private constructor, takes in an iterator over a vector of size_type
    IncidentIterator(std::vector<size_type>::iterator new_itr, 
		     Node root, const graph_type* graph_) : 
			incident_itr(new_itr), root_node(root), graph(graph_) {};
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
    EdgeIterator() {
    }

    /** Dereference operator for EdgeIterator
     * @pre this->edge_iter != edge_iter.end()
     * @return Edge with index @a i
    */
    Edge operator*() const {
        size_type index = edge_iter->first;
        return graph->edge(index);
    }
  
    /** Next operator for EdgeIterator
     * @return Reference to Edge Iterator with edge_iter pointing to next elem
    */
    EdgeIterator& operator++() {
	edge_iter++;
	return *this;
    }
   
    /** Equality operator for EdgeIterator
     * @return bool, whether Edge iterator has same graph, points to same edge
    */
    bool operator==(const EdgeIterator& other_iter) const {
	return (edge_iter == other_iter.edge_iter) && (graph == other_iter.graph);
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    std::map<size_type, std::pair<size_type, size_type>>::const_iterator edge_iter;
    const graph_type* graph; 

    /** Private constructor for edge iterator */
    EdgeIterator(std::map<size_type, std::pair<size_type, size_type>>::const_iterator iter,
		 const graph_type* graph_) :
			edge_iter(iter), graph(graph_) {} ;    
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns an iterator pointing to beginning of edge keys (set)*/
  edge_iterator edge_begin() const{
      return EdgeIterator(edge_keys.begin(), this);
  }
  
  /** Returns an iterator pointing to one past end of edge keys (set) */
  edge_iterator edge_end() const {
      return EdgeIterator(edge_keys.end(), this);
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    size_type uid;
    Point point; 
    node_value_type value;
    std::vector<size_type> connections;
  };

  std::map<size_type, internal_node> nodes;
  std::set<size_type> node_keys;

  //The pair containing the size_types will hold the indices
  //of the edges. uid might be redundant, but keeping it for now
  struct internal_edge {
    size_type uid; 
    Node a;
    Node b;
  }; 

  std::map<std::pair<size_type,size_type>, internal_edge> edges;
  std::map<size_type, std::pair<size_type, size_type>> edge_keys;

};

#endif // CME212_GRAPH_HPP
