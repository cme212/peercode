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

//template the Graph class
template <typename V, typename E>
class Graph {
 private:

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  //for templating of Graph class
  using node_value_type = V;
  using edge_value_type = E;

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

  //forward declare
  struct neighborNode;
  struct internalNode;
  struct internalEdge;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  //Graph data attributes
  size_type totnodes,totedges;
  std::vector<internalNode> nodevec;
  std::vector<internalEdge> edgevec;

  //unique ID vectors for removal on HW2
  std::vector<int> uniqNodeID, uniqEdgeID;

  Graph() {
    totnodes = 0;
    totedges = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //citing Graph_1038...see hw1 submission for details
  struct neighborNode {
      size_type edgeID;  //edge connecting the nodes
      size_type adjacentNode;
      //constructor for add_edge
      neighborNode(size_type id, size_type nextNode) {
          edgeID = id;
          adjacentNode = nextNode;
      }
  };


  //to carry all the 'heavy' items of a node and keep Node class light
  struct internalNode {
      Point nodepos;
      size_type nodeID;
      node_value_type value;
      //citing Graph_1038
      std::vector<neighborNode> neighborNodevec;
      //Constructor
      internalNode(Point p, size_type index, node_value_type val, std::vector<neighborNode> neighvec) {
          nodepos = p;
          nodeID = index;
          value = val;
          neighborNodevec = neighvec;
      }
  };

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
    }

    /** Return this node's position. */
    const Point& position() const {
      //geeksforgeeks help on -> used to access class members using pointers
      return mainGraph->nodevec[nodeindex].nodepos;
    }

    //Need a modifiable position method for HW2
    Point& position() {
      return mainGraph->nodevec[nodeindex].nodepos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return mainGraph->nodevec[nodeindex].nodeID;
    }

    // HW1: YOUR CODE HERE
    //return the value attribute from internalNode for that Node
    node_value_type& value() {
        return mainGraph->nodevec[nodeindex].value;
    }

    //return the value attribute from internalNode for that Node
    const node_value_type& value() const {
        return mainGraph->nodevec[nodeindex].value;
    }

    size_type degree() const {
        //get the size of the vector of neighbor nodes
        return mainGraph->nodevec[nodeindex].neighborNodevec.size();
    }

    incident_iterator edge_begin() const {
        //the center node is this node, the spoke starts with the 0 index
        //of the neighbor node vector
        return IncidentIterator(nodeindex,0,mainGraph);
    }

    incident_iterator edge_end() const {
        //need index of last 'spoke' node for the given node
        size_type lastspoke; //could maybe use degree function here?
        lastspoke = mainGraph->nodevec[nodeindex].neighborNodevec.size();
        return IncidentIterator(nodeindex,lastspoke,mainGraph);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //see if the graph and index are the same and return true
      if ((mainGraph == n.mainGraph && nodeindex == n.nodeindex)) {
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
      //make sure nodes are the same graph then check the index
      if (mainGraph<n.mainGraph || (mainGraph == n.mainGraph && nodeindex < n.nodeindex)) {
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* mainGraph;
    size_type nodeindex;
    //constructor to allow Graph class to construct node objects
    Node(const graph_type* graph, size_type indexID) {
        mainGraph = const_cast<graph_type*>(graph);
        nodeindex = indexID;
    }

  }; //last line of Node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //return nodevec.size(); //this was bad when adding remove
    return totnodes;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type nextvecindex;
    //get the next available index in the vector which is just its size
    nextvecindex = totnodes;
    //push_back can't accept everything that I need
    //stackoverflow used for emplace_back discussion accepting constructor
    //as arguments
    //citing Graph_1038 for HW1 additions to this method
    std::vector<neighborNode> theseAdj {};
    nodevec.emplace_back(position, nextvecindex,value,theseAdj);
    //HW2 for unique ID vector to perserve ID for removed nodes
    uniqNodeID.push_back(nodevec.size()-1);
    //the size has increased, so updates the number of nodes
    totnodes++;
    return Node(this, nodevec.size()-1); //didn't like mainGraph, 'this' works
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //check the node's graph is this graph instance
    //and check the index in within bounds
    if (this == n.mainGraph && n.nodeindex < nodevec.size()) {
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
    return Node(this, uniqNodeID[i]);
  }
  //end of public Graph class methods for nodes

  //'heavy' edge elements to keep Edge class light
  struct internalEdge {
      //Node A,B; //too heavy and slow fixed for hw2 on
      size_type node1ID, node2ID;
      size_type edgeID;
      edge_value_type value;
      //Constructor
      internalEdge(size_type a, size_type b, edge_value_type val,  size_type index) {
          node1ID = a;
          node2ID = b;
          value = val;
          edgeID = index;
      }
  };

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
      //access the internalEdge data attribute
      //return mainGraph->edgevec[edgeindex].A;
      return Node(mainGraph,n1ID);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return mainGraph->edgevec[edgeindex].B;
      return Node(mainGraph,n2ID);
    }

    //Values for edge on HW2
    //return the value attribute from internalEdge for that Edge
    edge_value_type& value() {
        return mainGraph->edgevec[edgeindex].value;
    }
    //return the value attribute from internalEdge for that Edge
    const edge_value_type& value() const {
        return mainGraph->edgevec[edgeindex].value;
    }

    //edge length method for HW2
    //https://www.cplusplus.com/reference/complex/norm/
    double length() const {
        double L;
        L = norm(node1().position() - node2().position());
        return L;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //check graph and then check nodes in BOTH directions
      if (mainGraph == e.mainGraph &&
              ((n1ID == e.n1ID && n2ID == e.n2ID) ||
               (n2ID == e.n1ID && n2ID == e.n1ID))) {
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
      //added condition to pass gtest hw2
      if (mainGraph<e.mainGraph || (mainGraph == e.mainGraph && edgeindex < e.edgeindex)) {
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* mainGraph;
    size_type edgeindex;
    size_type n1ID, n2ID;
    //constructor to allow Graph class to construct node objects
    Edge(const graph_type* graph, size_type indexID, size_type node1, size_type node2) {
        mainGraph = const_cast<graph_type*>(graph);
        edgeindex = indexID;
        n1ID = node1;
        n2ID = node2;
    }
  };  //last line of the edge class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    //return edgevec.size(); //had to change this as well...so many segfaults
    return totedges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, uniqEdgeID[i], edgevec[i].node1ID, edgevec[i].node2ID);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //take a look at each node attribute of each edge in edgevec
    //this would be more efficient with a pair or map...O(1) just a thought
    /*for (size_type i=0; i<edgevec.size(); i++) {
        if ((edgevec[i].A == a && edgevec[i].B == b) ||
            (edgevec[i].B == a && edgevec[i].A == b)) {
            return true;
        }
    }*/ //Fixed below for faster runtime on hw2...finally!

    //more efficient to just check through one nodes neighbor nodes
    for (size_type i=0; i<nodevec[a.nodeindex].neighborNodevec.size(); i++) {
        if (b.nodeindex == nodevec[a.nodeindex].neighborNodevec[i].adjacentNode) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
      size_type nextvectorindex;
      nextvectorindex = totedges;
      //use faster has_edge to check if the edge exists already
      if (has_edge(a,b) == false) {
          //then add a new edge
          edgevec.emplace_back(a.nodeindex, b.nodeindex, value, nextvectorindex);
          //add to the unique ID vector
          uniqEdgeID.push_back(nextvectorindex);
          //citing Graph_1038
          //for node a, add b as a connected node
          nodevec[a.nodeindex].neighborNodevec.push_back(neighborNode(\
	      edgevec.size()-1,b.nodeindex));
          //for node b, add a as a connected node
          nodevec[b.nodeindex].neighborNodevec.push_back(neighborNode(\
              edgevec.size()-1,a.nodeindex));
          totedges++;
      }
    return Edge(this, edgevec.size()-1, a.nodeindex, b.nodeindex);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //clear out each vector and set the counters to zero
    nodevec.clear();
    edgevec.clear();
    uniqNodeID.clear();
    uniqEdgeID.clear();
    totnodes = 0;
    totedges = 0;
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
    Node operator*() const {
        //need the index of the node being dereferenced using the iter index
        size_type indexref;
        indexref = nIterGraph->nodevec[nIterIndex].nodeID;
        return Node(nIterGraph, indexref);
    }
    NodeIterator& operator++() {
        size_type iterGraphsize;
        iterGraphsize = nIterGraph->num_nodes();
        //check that its not at the last element
        if (nIterIndex < iterGraphsize) {
            nIterIndex++;
        }
        //passed by reference and thus was modified, return the same iterator
        return *this;
    }
    bool operator==(const NodeIterator& nIter) const {
        if ((nIterGraph == nIter.nIterGraph) && (nIterIndex == nIter.nIterIndex)) {
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type nIterIndex;
    //didn't pass smoke test without const
    const graph_type* nIterGraph;

    //node iterator constructor
    NodeIterator(size_type index, const graph_type* graphptr) {
        nIterIndex = index;
        nIterGraph = graphptr;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
      //return the index 0 interator of this graph
      return NodeIterator(0,this);
  }
  node_iterator node_end() const {
      //return the index+1 which is just the size of the graph
      size_type lastelement;
      lastelement = this->num_nodes();
      return NodeIterator(lastelement,this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator>{
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
    Edge operator*() const {
        //return an edge object with the iter graph and both nodes
        size_type edgeID, node2;
        edgeID = iIterGraph->nodevec[centerNode].neighborNodevec[spokeNode].edgeID;
        node2 = iIterGraph->nodevec[centerNode].neighborNodevec[spokeNode].adjacentNode;
        return Edge(iIterGraph, edgeID, centerNode, node2);
    }

    IncidentIterator& operator++() {
        //move to the next neighbor node on the 'spoke'
        spokeNode++;
        return *this;
    }

    bool operator==(const IncidentIterator& iIter) const {
        if ((iIterGraph == iIter.iIterGraph) &&
            (centerNode == iIter.centerNode) &&
            (spokeNode == iIter.spokeNode)) {
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    //need graph pointer, node 1, and node 2 attributes
    size_type centerNode;
    size_type spokeNode;
    const graph_type* iIterGraph;
    //constructor
    IncidentIterator(size_type node1, size_type node2, const graph_type* graph) {
        centerNode = node1;
        spokeNode = node2;
        iIterGraph = graph;
    }
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
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
        //need the index of the edge being dereferenced using the iter index
        size_type node1, node2;
        //size_type indexref = eIterGraph->edgevec[eIterIndex].edgeID;
        node1 = eIterGraph->edgevec[eIterIndex].node1ID;
        node2 = eIterGraph->edgevec[eIterIndex].node2ID;
        return Edge(eIterGraph, eIterIndex, node1, node2);
    }

    EdgeIterator& operator++() {
            eIterIndex++;
        return *this;
    }

    bool operator==(const EdgeIterator& eIter) const {
        if ((eIterGraph == eIter.eIterGraph) && (eIterIndex == eIter.eIterIndex)) {
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type eIterIndex;
    const graph_type* eIterGraph;
    //constructor
    EdgeIterator(size_type index, const graph_type* graphptr) {
        eIterIndex = index;
        eIterGraph = graphptr;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
      return EdgeIterator(0,this);
  }
  edge_iterator edge_end() const {
      size_type lastelement;
      lastelement = this->num_edges();
      return EdgeIterator(lastelement,this);
  }

  //
  //HW2 remove methods section!
  //

  /** Remove an edge to the graph.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if an edge was removed, 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post new num_edges() == old num_edges() - 1.
   *
   * Complexity: No more than O(degree of node1^2 + num_edges())
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
      //first check if the edge even exists or maybe was removed
      if (has_edge(n1, n2) == false) {
          //post 377 from Ed...1 if an edge removed, 0 if not
          return 0;
      }
      //loop through the first nodes neighbors and look for the other node
      for (size_type i=0; i<n1.degree(); i++) {
          if (nodevec[n1.nodeindex].neighborNodevec[i].adjacentNode == n2.nodeindex) {
              //found our edge to delete, use erase() which requires an iter arg
              size_type eID = nodevec[n1.nodeindex].neighborNodevec[i].edgeID;
              nodevec[n1.nodeindex].neighborNodevec.erase(nodevec[n1.nodeindex].neighborNodevec.begin() + i);
              //now we need to remove the nodes from their neighbor node list
              for (size_type j=0; j<n2.degree(); j++) {
                  if (nodevec[n2.nodeindex].neighborNodevec[j].adjacentNode == n1.nodeindex) {
                      nodevec[n2.nodeindex].neighborNodevec.erase(nodevec[n2.nodeindex].neighborNodevec.begin() + j);
                  }
              }
              //erase the unique ID of the edge being deleted
              uniqEdgeID.erase(uniqEdgeID.begin() + edgevec[eID].edgeID);
              //update the edge id now for all edges after the one deleted
              for (size_type k=edgevec[eID].edgeID; k<uniqEdgeID.size(); k++) {
                  edgevec[uniqEdgeID[k]].edgeID = k;
              }
              //shrink the graph size by 1...iterator will still work
              totedges += -1;
          }
      }

      /* First attempt kept for completeness, will delete on hw3
      //get the node IDs
      size_type n1ID, n2ID;
      n1ID = n1.nodeindex;
      n2ID = n2.nodeindex;
      for (auto i=uniqEdgeID.begin(); i != uniqEdgeID.end(); ++i) {
          //find the corresponding edge
          if ((n1ID == edgevec[uniqEdgeID[*i]].node1ID &&
               n2ID == edgevec[uniqEdgeID[*i]].node2ID) ||
              (n1ID == edgevec[uniqEdgeID[*i]].node2ID &&
               n2ID == edgevec[uniqEdgeID[*i]].node1ID)) {
              //need to remove unique ID at the i index
              uniqEdgeID.erase(i);
              //need to decrement the tot number of edges
              totedges -= 1;
              //need to shift edge IDs after the removed edge by 1
              //j=i instead i+1 because we have already perfomed the erase
              for (size_type j = (*i); j<uniqEdgeID.size(); j++) {
                  edgevec[uniqEdgeID[j]].edgeID -= 1;
              }
              //need to update the neighbor node for each node
              for (size_type k=0; k != nodevec[n1ID].neighborNodevec.size(); k++) {
                  if (nodevec[n1ID].neighborNodevec[k].adjacentNode == n2ID) {
                      nodevec[n1ID].neighborNodevec.erase(nodevec[n1ID].neighborNodevec.begin() + k);
                  }
              }
              for (size_type k=0; k != nodevec[n2ID].neighborNodevec.size(); k++) {
                  if (nodevec[n2ID].neighborNodevec[k].adjacentNode == n1ID) {
                      nodevec[n2ID].neighborNodevec.erase(nodevec[n2ID].neighborNodevec.begin() + k);
                  }
              }
          }
          //return 1;
      }*/
      return 1;
  }

  //just use above in the following wrapper functions
  size_type remove_edge(const Edge& e) {
      node_type n1 = e.node1();
      node_type n2 = e.node2();
      return remove_edge(n1, n2);
  }
  edge_iterator remove_edge(edge_iterator e_it) {
      edge_type e = *e_it;
      //increment to the next before we erase
      //edge_iterator next_iter;
      //next_iter = (++e_it); //this idea doesn't work
      remove_edge(e);
      return e_it;
  }


  /** Remove a node from the graph and delete any connected edges.
   * @pre @a a is a distinct valid nodes of this graph
   * @return a 1 if node removed, 0 otherwise.
   * @post has_node(@a a) == false
   * @post new num_nodes() == old num_nodes() - 1.
   *
   * Complexity: No more than O(degree of node + remove_edge complexity)
   */
  size_type remove_node(const Node& n) {
      if (has_node(n) == false) {
          return 0;
      }
      size_type nID = n.index();  //n.nodeindex;
          //need to remove any edge spawning from that node
          //for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
          //    remove_edge(*iter);
          //}
          //this attempt continued to segfault, some peers recommened a while
          //loop with the degree function essentially counting down as edges
          //are removed
          while (n.degree() > 0) {
              auto deledge = *(n.edge_begin());
              remove_edge(deledge);
          }
          //need to remove that node from the unique ID list
          //erase function needs an iterator not an index
          uniqNodeID.erase(uniqNodeID.begin() + nID);
          //need to decrement the graph size by 1
          totnodes -= 1;
          //need to decrement every node index after by 1
          for (size_type i = nID; i<uniqNodeID.size(); i++) {
              nodevec[uniqNodeID[i]].nodeID = i; //-= 1;
          }
      return 1;
  }

  node_iterator remove_node(node_iterator n_it) {
      node_type n = *n_it;
      //increment to the next before we erase
      //node_iterator next_iter;
      //next_iter = (++n_it);
      remove_node(n);
      return n_it;
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
