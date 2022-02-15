#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * 
 * NOTE: Inspiration taken from graph_976.hpp peercode, 2022.10.26
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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

  using node_value_type = V;
  typedef E edge_value_type; // Different form of notation, means the same thing

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
  using length_type = double;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    nid_node_map = {};
    eid_edge_map = {};
    edge_eid_map = {};
    nid_nidx = {};
    eid_eidx = {};
    nidx_nid = {};
    eidx_eid = {};
    nid_count = 0;
    eid_count = 0;
  }

  /** Default destructor */
  ~Graph() = default;















  //
  // NODES
  //
  // class Node : private totally_ordered<Node>;

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>{
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

    // Node(Graph* this_graph, size_type this_nid) {
    //   // instantiate new node with pointer to graph and index
    //   my_graph = this_graph;
    //   nid = this_nid;
    // }

    Node(const Graph* this_graph, size_type this_nid) 
    : my_graph(const_cast<Graph*>(this_graph)), nid(this_nid) {
    }

    //     SimpleElement(const Ex_Graph* set, size_type uid)
    //     : set_(const_cast<Ex_Graph*>(set)), uid_(uid) {
    // }

    Node() {}

    /** Return this node's value. */
    node_value_type& value() {
      // get value
      return this->my_graph->nid_node_map.at(this->nid).my_value;
    }
    
    /** Return this node's value. */
    const node_value_type& value() const {
      // get value
      return this->my_graph->nid_node_map.at(this->nid).my_value;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // get internal node
      return this->my_graph->nid_node_map.at(this->nid).my_point;
    }

    /** Return this node's position. */
    Point& position() {
      // HW0: YOUR CODE HERE
      // get internal node
      return this->my_graph->nid_node_map.at(this->nid).my_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->my_graph->nid_nidx[this->nid];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Returns the number of edges attached to the node */
    size_type degree() const {
      // Get size of unordered map in 
      if (this->my_graph->edge_eid_map.count(this->nid) == 1) {
        return this->my_graph->edge_eid_map.at(this->nid).size();
      } else {
      return 0;
      }
    }

    /** Returns the first edge connected to this node 
     * Note: Underlying datastructure is a map, hence the verbosity
     * Primarily relies on built in map iterator
    */
    incident_iterator edge_begin() const {
      IncidentIterator incident_it = IncidentIterator();
      incident_it.my_graph = this->my_graph;
      incident_it.my_nid = this->nid;
      // In case where no connecting nodes, make dummy map_iter
      if (Node(my_graph, nid).degree() == 0) {
        std::unordered_map<size_type, size_type> dummy = {};
        incident_it.map_iter = dummy.begin();
        // Otherwise, use actual map
      } else {
        incident_it.map_iter = this->my_graph->edge_eid_map.at(nid).begin();
      }
      return incident_it;
    }

    /** Returns one past the final connected to this node 
     * Note: Underlying datastructure is a map, hence the verbosity
    */
    incident_iterator edge_end() const {
      IncidentIterator incident_it = IncidentIterator();
      incident_it.my_graph = this->my_graph;
      incident_it.my_nid = this->nid;
      // In case where no connecting nodes, make dummy map_iter
      if (Node(my_graph, nid).degree() == 0) {
        std::unordered_map<size_type, size_type> dummy = {};
        incident_it.map_iter = dummy.end();
      } else {
        incident_it.map_iter = this->my_graph->edge_eid_map.at(nid).end();
      }
      return incident_it;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool graph_equal = (this->my_graph == n.my_graph);
      bool nid_equal = (this->nid == n.nid);       
      return (graph_equal & nid_equal);
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
      // Compare indices, see if equal
      if (this->my_graph == n.my_graph) {
        return this->nid < n.nid;
      } else {
        // otherwise, return a comparison of the graph addresses
        return std::addressof(this->my_graph) < std::addressof(n.my_graph);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // the node has a pointer to graph 
    Graph* my_graph;
    // the node has an index
    size_type nid;
  };



  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    size_type my_size = this->nid_node_map.size();
    return my_size;
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
  // called on the graph object, self
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
      // HW0: YOUR CODE HERE
      // Create internal node, add to set 
      size_type nidx = num_nodes();
      size_type nid = nid_count;
      // Increment unique node id count
      nid_count = nid_count + 1;
      // (in order to not copy over pointer)
      internal_node my_int_node;
      my_int_node.my_point = Point(position.x, position.y, position.z);
      // assign value to my_value in internal_node
      my_int_node.my_value = value;
      // add to nid -> internal_node map
      nid_node_map[nid] = my_int_node;
      // add to nid_nidx map and nidx_nid map
      nid_nidx[nid] = nidx;
      nidx_nid.push_back(nid);
      // Create node and return
      Node my_node = Node(this, nid);
      return my_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // same nid and graph is sufficient for containing node 
    size_type node_id = n.nid;
    if(nid_nidx.find(node_id) != nid_nidx.end()) {
      if (n.my_graph == this) {
        return true;
      }
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
    assert(i < num_nodes());
    return Node(this, nidx_nid.at(i));        
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] _a_ the node to be removed, must be in this graph
   * @return 1 if sucessful, 0 if not successful
   * @post has_node(_a_) == false
   * @post new num_nodes() == old num_nodes() - 1 if sucessful, otherwise no change
   * @post all edges connected to the node are removed as well
   * 
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   * Will remove edges connected to node, old edge with index j may not 
   * equal new edge at index j.  Some edges may no longer be valid.
   * Removing node may invalidate outstanding edge iterators (including incident
   * iterator) if they point at an edge that has been removed.
   * Node iterators may also be invalidated
   * 
   *
   * Complexity: O(1), in average case
   * To achive constant runtime, we first move the edge to be removed to the end
   * of the nid_node storage vector (by swapping with element currently at end),
   * and then pop last element.
   */
  size_type remove_node(const Node& a) {
    // Check if node is from same graph 
    if (this != a.my_graph) {
      return 0;
    }
    // Check that node in graph
    if (!has_node(a)) {
      return 0;
    }
    // Remove any incident edges
    // CONCERN: check for infinite loop!
    while (a.edge_begin() != a.edge_end()) {
      // Remove the first edge 
      remove_edge(*a.edge_begin());
    }
    // Save important information
    size_type remove_nid = a.nid;
    size_type remove_nidx = nid_nidx.at(remove_nid);
    size_type swap_nidx = num_nodes() -1;
    size_type swap_nid = nidx_nid.at(swap_nidx);
    
    // Remove from datastructures
    // (1) nid_node_map
    nid_node_map.erase(remove_nid);
    // (2) Perform swap, removing from nid_nidx unordered map, nidx_nid vector
    nid_nidx.at(swap_nid) = remove_nidx;
    nid_nidx.erase(remove_nid);
    nidx_nid.at(remove_nidx) = swap_nid;
    nidx_nid.pop_back();
    // (3) ensure nid_count unchanged

    return 1;
  }


   /** Add a node to the graph, returning the added node.
   * @param[in] @ n_it iterator to the node to be removed, must be from this graph
   * @return iterator pointing to element that would follow this one, n_it++
   * @post has_node(*n_it)) == false
   * @post new num_nodes() == old num_nodes() - 1 if sucessful, otherwise no change
   * @post all edges connected to the node are removed as well
   * 
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   * Will remove edges connected to node, old edge with index j may not 
   * equal new edge at index j.  Some edges may no longer be valid.
   * Removing node may invalidate outstanding edge iterators (including incident
   * iterator) if they point at an edge that has been removed.
   * Other node iterators may also be invalidated.
   * 
   * Complexity: O(1), in average case
   */
  node_iterator remove_node(node_iterator n_it) {
    // Call to remove_node on the node
    size_type success = remove_node(n_it);
    // If successful, we removed a node and the current iterator points 
    // to the next node.  Else, we return incremented node.
    if (success) {
      return n_it;
    } else {
      return n_it++;
    }
  }





















  //
  // EDGES  legooooo
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:

    /** Construct a valid edge to return **/
    Edge(const Graph* this_graph, size_type this_eid) 
    : my_graph(const_cast<Graph*>(this_graph)), eid(this_eid) {
    }

    /** Construct an invalid Edge. */
    Edge() {
      // See if deleting this is okay!
    }

    /** Return reference to value stored in this edge */
    edge_value_type& value() {
      // go to graph to get interal edge and retrieve value
      return this->my_graph->eid_edge_map.at(this->eid).my_value;
    }

     /** Return reference to value stored in this edge */
    const edge_value_type& value() const {
      // go to graph to get interal edge and retrieve value
      return this->my_graph->eid_edge_map.at(this->eid).my_value();
    }

    /** Return the length of the edge, distance between connected nodes */
    length_type length() {
      // return current distance between two nodes
      return norm(this->node1().position() - this->node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return this->my_graph->eid_edge_map.at(this->eid).node_1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return this->my_graph->eid_edge_map.at(this->eid).node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE      
      bool same_graph = (this->my_graph == e.my_graph);
      // case 1: n1=n1, n2=n2
      bool node_same = ((this->node1() == e.node1()) & \
      (this->node2() == e.node2()));
      bool node_flipped = ((this->node1() == e.node2()) & \
      (this->node2() == e.node1()));
      
      return (same_graph & (node_flipped | node_same));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // if from same graph, compare unique ids of the edges
      if (this->my_graph == e.my_graph) {
        return this->eid < e.eid;
      } else {
        // otherwise, return a comparison of the graph addresses
        return std::addressof(this->my_graph) < std::addressof(e.my_graph);
      }
    }

    /** Return the index of this edge 
     * 
     */
     size_type index() const {
         return this->my_graph->eid_eidx[this->eid];
     } 

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* my_graph;
    size_type eid;
  };



  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->eid_edge_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // assert we are asking for valid edge
    assert(i < num_edges());
    // return new node
    return Edge(this, eidx_eid.at(i));        
  }

  

 /** Returns the index of the edge connecting two valid nodes.
  * @pre assumes valid distinct nodes
  * @return index if valid, 4200000000 if not valid **/
  size_type index_eid(const Node& a, const Node& b) const {
      // HW0: my code
      // Any ordering works!
      if (a.my_graph->edge_eid_map.count(a.nid) == 1) {
          if (a.my_graph->edge_eid_map.at(a.nid).count(b.nid) == 1) {
              return a.my_graph->edge_eid_map.at(a.nid).at(b.nid);
          }
      }
    return 4200000000;
  }

/** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // see if index is invalid
    if (index_eid(a,b) == 4200000000) {
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
    assert(a != b);
    // HW0: YOUR CODE HERE
    // check if edge already exists 
    if (has_edge(a,b)) {
      // Make a new edge
      internal_edge my_int_edge;
      // add to eid_edge_map
      my_int_edge.node_1 = a;
      my_int_edge.node_2 = b;
      a.my_graph->eid_edge_map[index_eid(a,b)] = my_int_edge;
    // return the current edge
    return Edge(this,index_eid(a,b));
    };
    // decleare index and new internal edge
    size_type eidx = num_edges();
    size_type eid = eid_count;
    eid_count = eid_count + 1;
    // Add to eid and eidx conversion datastructures
    eid_eidx[eid] = eidx;
    eidx_eid.push_back(eid);
    internal_edge my_int_edge;
    // add to eid_edge_map
    my_int_edge.node_1 = a;
    my_int_edge.node_2 = b;
    a.my_graph->eid_edge_map[eid] = my_int_edge;
    // add to edge_eid_map
    // if there's not already an unordered map, make one
    if (a.my_graph->edge_eid_map.count(a.nid) == 0) {
        a.my_graph->edge_eid_map[a.nid] = {};
    }
    a.my_graph->edge_eid_map.at(a.nid)[b.nid] = eid;

    // if there's not already an unordered map, make one
    if (b.my_graph->edge_eid_map.count(b.nid) == 0) {
        b.my_graph->edge_eid_map[b.nid] = {};// 
    }
    b.my_graph->edge_eid_map.at(b.nid)[a.nid] = eid;

    // return the new edge
    return Edge(this, eid);        
  }

  /** Remove an edge from the graph.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if sucessfully removed, 0 if node not present (and thus not removed)
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). May invalidate edge iterators.
   *
   * Complexity: O(1), in average case
   * To achive constant runtime, we first move the edge to be removed to the end
   * of the eid_edge storage vector (by swapping with element currently at end),
   * and then pop last element.
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // If node from another graph, return 0 signaling fialure
    if( (a.my_graph != this) | (b.my_graph != this)) {
      return 0;
    }
    // Get unique edge id (aka eid)
    size_type remove_eid = index_eid(a,b);
    // Check if in graph
    if (remove_eid == 4200000000) {
      return 0;
    }
    // Save key values 
    size_type remove_eidx = eid_eidx.at(remove_eid);
    size_type swap_eidx = num_edges() - 1;
    size_type swap_eid = eidx_eid.at(swap_eidx);
    
    // Remove from all datastructures
    // (1) edge_eid_map (maps edge by nodes to edge id)
    edge_eid_map.at(a.nid).erase(b.nid);
    edge_eid_map.at(b.nid).erase(a.nid);
    // (2) eid_edge_map (stores internal_edge by eid)
    eid_edge_map.erase(remove_eid);
    // (3) eid_eidx (unordered map) and eidx_eid (vector)
    eid_eidx.at(swap_eid) = remove_eidx; //map swap edge eid to eidx of removed node
    eid_eidx.erase(remove_eid); // remove eid of removed edge  
    eidx_eid.at(remove_eidx) = swap_eid;
    eidx_eid.pop_back(); // remove last entry in eidx_eid vector
    // (4) note eid_count is unchanged
    
    return 1;
  }

  /** Remove an edge from the graph.
   * @pre _e_ is a valid edge in graph, with nodes _a_ and _b_
   * @return 1 if sucessfully removed, 0 if node not present (and thus not removed)
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). May invalidate edge iterators.
   *
   * Complexity: O(1), in average case
   * To achive constant runtime, we first move the edge to be removed to the end
   * of the eid_edge storage vector (by swapping with element currently at end),
   * and then pop last element.
   */
  size_type remove_edge(const Edge& e) {
    // If from another graph, return false
    if (e.my_graph != this){
      return 0;
    }
    // check to see if not in graph 
    if (eid_edge_map.count(e.eid) == 0) {
      return 0;
    }
    // Call remove_edge on nodes 
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph, based on iterator.
   * @pre _e_it_ is a valid edge iterator for edge between nodes _a_ and _b_
   * @return iterator pointing to next edge to evaluate (already incremetned)
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). May invalidate edge iterators.
   *
   * Complexity: O(1), in average case
   * To achive constant runtime, we first move the edge to be removed to the end
   * of the eid_edge storage vector (by swapping with element currently at end),
   * and then pop last element.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // assert the graph of e_it is the same as this one
    assert(e_it.my_graph == this);
    // get edge object from edge iterator 
    edge_type e = *e_it;
    // remove edge
    remove_edge(e);
    // return iterator pointing at next edge, which is the same iterator
    // since the edge iterator is based on index of elements 
    // (can also design implementation to return edge_begin()) but this 
    // implementation choice facilitates the ability to iterate over nodes,
    // and remove with some consistency
    return e_it;
  }




  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // clear all datastructures
    this->nid_node_map = {};
    this->eid_edge_map = {};
    this->edge_eid_map = {};
    this->nid_nidx = {};
    this->nidx_nid = {};
    this->eid_eidx = {};
    this->eidx_eid = {};
    this->nid_count = 0;
    this->eid_count = 0;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct a valid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    //Increments to the next node in node class.
    NodeIterator& operator++() {
      this->my_nidx += 1;
      return *this;
    }
    
    /** Defines equality between two iterators */
    bool operator==(const NodeIterator& node_iter) const {
        // check if my_graph, my_index is the same
        return ((this->my_graph == node_iter.my_graph) & \
        (this->my_nidx == node_iter.my_nidx));      
        }
    
    /** Defines inequality between two iterators */
    bool operator!=(const NodeIterator& node_iter) {
        // check if my_graph, my_index is the same
        return !(*this == node_iter);      
        }

    /** Dereferences iterator */
    Node operator*() const {
        // Returns the Node at this index
        return Node(this->my_graph,this->my_graph->nidx_nid.at(this->my_nidx));
    }

   private:
    friend class Graph;
    // Store the graph and index of nodes
    const graph_type* my_graph;
    size_type my_nidx;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /** Make beginning of node iterator */
    NodeIterator node_begin() const {
      NodeIterator it = NodeIterator();
      it.my_graph = this;
      it.my_nidx = 0;
      return it;
    }

    /** Make end of node iterator */
    NodeIterator node_end() const {
      NodeIterator it = NodeIterator();
      it.my_graph = this;
      it.my_nidx = num_nodes(); // one past the number of nodes
      return it;
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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereference the iterator */
    Edge operator*() const{
      size_type other_nid = map_iter->first;
      size_type eid = map_iter->second;
      internal_edge this_int_edge = my_graph->eid_edge_map.at(eid);
      // MAKE SURE MODIFY EDGE TO HAVE node1 = my_node, node2 = other_node
      if (this_int_edge.node_1.nid != my_nid) {
        // Swap nodes in this_int_edge
        this_int_edge.node_1 = Node(my_graph, my_nid);
        this_int_edge.node_2 = Node(my_graph,other_nid);
      }
      // MUST REASSIGN INTERNAL_EDGE TO EID_EDGE_MAP (!!)
      my_graph->eid_edge_map[eid] = this_int_edge;
      return Edge(my_graph, eid);
    }

    /** Increment the iterator */
    IncidentIterator& operator++() {
      // Increment mapiter
      map_iter++;
      return *this;
    }

    /** Compare to see if equal 
     * Note: != method implemented automatically by totally_ordered
     * */
    bool operator==(const IncidentIterator& other_it) const {
      // compare graphs
      bool same_graph = (this->my_graph == other_it.my_graph);
      // compare nodes
      bool same_nodes = (this->my_nid == other_it.my_nid);
      // compare node iterators
      bool same_iter = (this->map_iter == other_it.map_iter);
      // return true if all true 
      return (same_graph & same_nodes & same_iter);
    }

   private:
    friend class Graph;
    graph_type* my_graph; // This graph
    size_type my_nid;  // Unique node id of node whose edges we are iterating over
    // C++ built in map iterator
    std::unordered_map<size_type, size_type>::iterator map_iter;

    // HW1 #3: YOUR CODE HERE

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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Dereference the iterator */
    Edge operator*() const {
      size_type my_eid = my_graph->eidx_eid.at(my_eidx);
      return Edge(my_graph, my_eid);
    }

    /** Increment the iterator */
    EdgeIterator& operator++() {
      my_eidx++;
      return *this;
    }

    /** Test equality of the iterator */
    bool operator==(const EdgeIterator& other_edge_it) const {
      bool same_graph = (this->my_graph == other_edge_it.my_graph);
      bool same_eidx = (this->my_eidx == other_edge_it.my_eidx);
      return (same_graph & same_eidx);
    }

    /** Test inequality of the iterator */
    bool operator!=(const EdgeIterator& other_edge_it) const {
      return !(*this == other_edge_it);
    }

   private:
    friend class Graph;
    // Iterator based on index of edges
    // Since Edge(i) has O(1) lookup time
    const graph_type* my_graph;
    size_type my_eidx;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns iterator pointing to first edge */
  edge_iterator edge_begin() const {
    EdgeIterator my_edge_it = EdgeIterator();
    my_edge_it.my_graph = this;
    my_edge_it.my_eidx = 0;
    return my_edge_it;
  }

  /** Returns iterator pointing to one past last edge */
  edge_iterator edge_end() const {
    EdgeIterator my_edge_it = EdgeIterator();
    my_edge_it.my_graph = this;
    my_edge_it.my_eidx = num_edges(); // One past the last 
    return my_edge_it;
  }







 private:

  /*  NODE DATASTRUCTURES */
  // internal_node DATASTRUCTURE
  // Contains the point corresponding to each node
  struct internal_node {
    Point my_point;
    node_value_type my_value;
  };

  // idx_node_map
  // Maps unique node id (nid) to internal_node
  std::unordered_map<size_type, internal_node> nid_node_map;

  // nid_nidx
  // maps from unique node id to node index, map 
  std::unordered_map<size_type,size_type> nid_nidx;

  // nidx_nid
  // 1D vector with elements of nid's in order
  std::vector<size_type> nidx_nid;

  // nid_count
  // keep track of next nid to use, to ensure unique
  size_type nid_count;





  /* EDGE DATASTRUCTURES */
  // internal_edge DATASTRUCTURE
  // Contains both nodes
    struct internal_edge {
      Node node_1;
      Node node_2;
      edge_value_type my_value;
  };

  // eid_edge_map
  // maps unique edge id (eid) to internal_edge datastructure
  std::unordered_map<size_type,internal_edge> eid_edge_map;

  // edge_eid_map
  // maps node_1 to map of connected nodes, with corresponding eids
  std::unordered_map<size_type, \
  std::unordered_map<size_type, size_type>> edge_eid_map;

  // eid_eidx
  // maps from unique edge id to edge index, map 
  std::unordered_map<size_type,size_type> eid_eidx;

  // eidx_eid
  // 1D vector with elements of eid's in order
  std::vector<size_type> eidx_eid;

  // eid_count
  // keep track of next eid to use, to ensure unique
  size_type eid_count;

};

#endif // CME212_GRAPH_HPP